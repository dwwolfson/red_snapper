
library(shiny)
library(here)
library(dplyr)
library(broom)
library(dotwhisker)


# call all the processed data
glm_5m <- readRDS(here("processed_data","glm_ready_3fish_5_m^2_res.Rdata"))
glm_10m <- readRDS(here("processed_data","glm_ready_3fish_10_m^2_res.Rdata"))
glm_15m <- readRDS(here("processed_data","glm_ready_3fish_15_m^2_res.Rdata"))
glm_20m <- readRDS(here("processed_data","glm_ready_3fish_20_m^2_res.Rdata"))
glm_25m <- readRDS(here("processed_data","glm_ready_3fish_25_m^2_res.Rdata"))
glm_30m <- readRDS(here("processed_data","glm_ready_3fish_30_m^2_res.Rdata"))
glm_35m <- readRDS(here("processed_data","glm_ready_3fish_35_m^2_res.Rdata"))
glm_40m <- readRDS(here("processed_data","glm_ready_3fish_40_m^2_res.Rdata"))
glm_45m <- readRDS(here("processed_data","glm_ready_3fish_45_m^2_res.Rdata"))
glm_50m <- readRDS(here("processed_data","glm_ready_3fish_50_m^2_res.Rdata"))

# fit models for each ID
by_snap_5m <- glm_5m %>%
    group_by(id) %>%                        # group data by transmission
    do(tidy(glm(z~habitat+crw, offset=log(tau), family=poisson, data=.))) %>%         # run model on each group
    ungroup() %>% rename(model = id) %>% mutate(size = 5)

by_snap_10m <- glm_10m %>%
    group_by(id) %>%                        # group data by transmission
    do(tidy(glm(z~habitat+crw, offset=log(tau), family=poisson, data=.))) %>%         # run model on each group
    ungroup() %>% rename(model = id) %>% mutate(size = 10)

by_snap_15m <- glm_15m %>%
    group_by(id) %>%                        # group data by transmission
    do(tidy(glm(z~habitat+crw, offset=log(tau), family=poisson, data=.))) %>%         # run model on each group
    ungroup() %>% rename(model = id) %>% mutate(size = 15)

by_snap_20m <- glm_20m %>%
    group_by(id) %>%                        # group data by transmission
    do(tidy(glm(z~habitat+crw, offset=log(tau), family=poisson, data=.))) %>%         # run model on each group
    ungroup() %>% rename(model = id) %>% mutate(size = 20)

by_snap_25m <- glm_25m %>%
    group_by(id) %>%                        # group data by transmission
    do(tidy(glm(z~habitat+crw, offset=log(tau), family=poisson, data=.))) %>%         # run model on each group
    ungroup() %>% rename(model = id) %>% mutate(size = 25)

by_snap_30m <- glm_30m %>%
    group_by(id) %>%                        # group data by transmission
    do(tidy(glm(z~habitat+crw, offset=log(tau), family=poisson, data=.))) %>%         # run model on each group
    ungroup() %>% rename(model = id) %>% mutate(size = 30)

by_snap_35m <- glm_35m %>%
    group_by(id) %>%                        # group data by transmission
    do(tidy(glm(z~habitat+crw, offset=log(tau), family=poisson, data=.))) %>%         # run model on each group
    ungroup() %>% rename(model = id) %>% mutate(size = 35)

by_snap_40m <- glm_40m %>%
    group_by(id) %>%                        # group data by transmission
    do(tidy(glm(z~habitat+crw, offset=log(tau), family=poisson, data=.))) %>%         # run model on each group
    ungroup() %>% rename(model = id) %>% mutate(size = 40)

by_snap_45m <- glm_45m %>%
    group_by(id) %>%                        # group data by transmission
    do(tidy(glm(z~habitat+crw, offset=log(tau), family=poisson, data=.))) %>%         # run model on each group
    ungroup() %>% rename(model = id) %>% mutate(size = 45)

by_snap_50m <- glm_50m %>%
    group_by(id) %>%                        # group data by transmission
    do(tidy(glm(z~habitat+crw, offset=log(tau), family=poisson, data=.))) %>%         # run model on each group
    ungroup() %>% rename(model = id) %>% mutate(size = 50)

# rbind all the grid size coefficient values
by_snap <- rbind(by_snap_5m, by_snap_10m, by_snap_15m, by_snap_20m, by_snap_25m, by_snap_30m, by_snap_35m, by_snap_40m, by_snap_45m, by_snap_50m)

by_snap$size = as.factor(by_snap$size)


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("CTMC GLM Coefficient values across grid sizes"),

    # Sidebar with a slider input for grid sizes
    sidebarLayout(
        sidebarPanel(
            sliderInput("size",
                        "Grid size (m)",
                        min = 5,
                        max = 50,
                        value = 5,
                        step = 5)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("Plot")
        )
    )
)

# Define a server for the shiny app
server <- function(input, output) {
    
    # Create a reactive dataset for plotting by grid size
    plotInput <- reactive({
        
        tmp <- by_snap[by_snap$size==input$"size",]
        
        small_multiple(tmp) +
            theme_bw() + ylab("Coefficient Estimate") +
            geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
            ggtitle("CTMC GLM") +
            theme(plot.title = element_text(face = "bold"), 
                  legend.position = "none",
                  axis.text.x = element_text(angle = 60, hjust = 1)) + xlab("Snapper ID") +
            scale_y_continuous(breaks = seq(-2,0.5,0.2)) + 
            coord_cartesian(ylim = c(-2, 0.5))
            })
    
    
    # Fill in the spot we created for a plot
    output$Plot <- renderPlot({
        plotInput()}, width =800, height = 550)
}

# Run the application 
shinyApp(ui = ui, server = server)
