library(shiny)
library(data.table)
library(bslib)
library(DT)
library(ggplot2)

ui <- page_sidebar(
  title = "RData File Upload",
  fillable = TRUE,  # Make the page fillable
  
  # Sidebar with all controls
  sidebar = sidebar(
    fileInput("file", "Upload .RData file",
              accept = c(".RData", ".rda")),
    helpText("Please upload an .RData file containing a data.table object"),
    uiOutput("object_selector"),
    uiOutput("x_var_selector"),
    uiOutput("y_var_selector"),
    uiOutput("color_var_selector"),
    uiOutput("threshold_slider")
  ),
  
  # Main panel with cards in a scrollable container
  div(
    style = "height: 100vh; overflow-y: auto;",  # Make container scrollable
    card(
      card_header("Scatter Plot"),
      plotOutput("scatter_plot", height = "500px")
    ),
    
    card(
      card_header("Color Variable Distribution"),
      plotOutput("color_histogram", height = "300px")
    ),
    
    card(
      card_header("Data Preview"),
      DT::dataTableOutput("table")
    ),
    
    card(
      card_header("Data Information"),
      textOutput("data_info"),
      verbatimTextOutput("structure")
    )
  )
)

server <- function(input, output, session) {
  
  # Reactive value to store loaded objects
  loaded_objects <- reactiveVal(NULL)
  selected_dt <- reactiveVal(NULL)
  
  # Handle file upload
  observeEvent(input$file, {
    req(input$file)
    
    # Create new environment to load data
    env <- new.env()
    
    tryCatch({
      load(input$file$datapath, envir = env)
      
      # Find all data.table objects in the environment
      dt_objects <- names(eapply(env, function(x) is.data.table(x)))
      
      if (length(dt_objects) == 0) {
        showNotification("No data.table objects found in the file", type = "error")
        return(NULL)
      }
      
      loaded_objects(env)
      
      # If only one data.table, select it automatically
      if (length(dt_objects) == 1) {
        selected_dt(get(dt_objects[1], envir = env))
      }
      
    }, error = function(e) {
      showNotification(paste("Error loading file:", e$message), type = "error")
    })
  })
  
  # Create dynamic object selector
  output$object_selector <- renderUI({
    req(loaded_objects())
    
    dt_objects <- names(eapply(loaded_objects(), function(x) is.data.table(x)))
    
    if (length(dt_objects) > 1) {
      selectInput("selected_object", "Choose data.table:", 
                  choices = dt_objects)
    }
  })
  
  # Create X variable selector
  output$x_var_selector <- renderUI({
    req(selected_dt())
    var_choices <- names(selected_dt())
    selectInput("x_var", "Select X Variable:",
                choices = var_choices,
                selected = var_choices[1])
  })
  
  # Create Y variable selector
  output$y_var_selector <- renderUI({
    req(selected_dt())
    var_choices <- names(selected_dt())
    selectInput("y_var", "Select Y Variable:",
                choices = var_choices,
                selected = var_choices[min(2, length(var_choices))])
  })
  
  # Create Color variable selector
  output$color_var_selector <- renderUI({
    req(selected_dt())
    var_choices <- names(selected_dt())
    selectInput("color_var", "Select Color Variable:",
                choices = var_choices,
                selected = var_choices[min(3, length(var_choices))])
  })
  
  # Create threshold slider
  output$threshold_slider <- renderUI({
    req(selected_dt(), input$color_var)
    color_values <- selected_dt()[[input$color_var]]
    if(is.numeric(color_values)) {
      sliderInput("threshold", "Color Variable Threshold:",
                  min = min(color_values, na.rm = TRUE),
                  max = max(color_values, na.rm = TRUE),
                  value = median(color_values, na.rm = TRUE))
    }
  })
  
  # Update selected data.table when user changes selection
  observeEvent(input$selected_object, {
    req(input$selected_object, loaded_objects())
    selected_dt(get(input$selected_object, envir = loaded_objects()))
  })
  
  # Create scatter plot
  output$scatter_plot <- renderPlot({
    req(selected_dt(), input$x_var, input$y_var, input$color_var)
    
    dt <- selected_dt()
    
    # Create factor based on threshold if numeric
    if(is.numeric(dt[[input$color_var]])) {
      req(input$threshold)
      color_factor <- factor(dt[[input$color_var]] > input$threshold,
                           labels = c("Below threshold", "Above threshold"))
    } else {
      color_factor <- factor(dt[[input$color_var]])
    }
    
    # Create transformed X and Y variables
    # x_transformed <- log(dt[[input$x_var]] - min(dt[[input$x_var]], na.rm = TRUE) + 1)
    # y_transformed <- log(dt[[input$y_var]] - min(dt[[input$y_var]], na.rm = TRUE) + 1)

    x_transformed <- sqrt(dt[[input$x_var]] - min(c(dt[[input$x_var]],
                                                   dt[[input$y_var]]),
                                                  na.rm = TRUE) + 1)
    y_transformed <- sqrt(dt[[input$y_var]] - min(c(dt[[input$x_var]],
                                                   dt[[input$y_var]]),
                                                  na.rm = TRUE) + 1)
    
    # Create temporary data frame with transformed variables
    plot_dt <- data.frame(
      x = x_transformed,
      y = y_transformed,
      color = color_factor
    )
    
    ggplot(plot_dt, aes(x = x, y = y)) +
      geom_point(aes(color = color), alpha = 0.6) +  # Changed alpha to 0.6
      theme_bw() +
      labs(x = paste("log(", input$x_var, " - min(", input$x_var, ") + 1)"),
           y = paste("log(", input$y_var, " - min(", input$y_var, ") + 1)"),
           color = input$color_var) +
      scale_color_manual(values = c("blue", "red")) +
      theme(legend.position = "bottom")
  })
  
  # Create histogram for color variable
  output$color_histogram <- renderPlot({
    req(selected_dt(), input$color_var)
    
    dt <- selected_dt()
    color_var <- dt[[input$color_var]]
    
    p <- ggplot(dt, aes_string(x = input$color_var)) +
      theme_bw() +
      labs(x = input$color_var, y = "Count")
    
    if(is.numeric(color_var)) {
      req(input$threshold)
      p <- p + 
        geom_histogram(aes(fill = color_var > input$threshold), 
                      bins = 30, alpha = 0.6) +  # Changed alpha to 0.6
        geom_vline(xintercept = input$threshold, 
                  color = "red", linetype = "dashed", size = 1) +
        scale_fill_manual(values = c("blue", "red"),
                         labels = c("Below threshold", "Above threshold"),
                         name = "Groups") +
        theme(legend.position = "bottom")
    } else {
      p <- p + 
        geom_bar(fill = "skyblue", alpha = 0.6) +  # Changed alpha to 0.6
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
    
    p
  })
  
  # Display the data.table
  output$table <- DT::renderDataTable({
    req(selected_dt())
    selected_dt()
  }, options = list(scrollX = TRUE))
  
  # Display information about the data
  output$data_info <- renderText({
    req(selected_dt())
    dt <- selected_dt()
    paste("Dimensions:", nrow(dt), "rows and", ncol(dt), "columns")
  })
  
  # Display structure of the data
  output$structure <- renderPrint({
    req(selected_dt())
    str(selected_dt())
  })
}

shinyApp(ui, server)
