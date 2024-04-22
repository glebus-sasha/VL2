library(shiny)
library(plotly)
library(dplyr)

# UI
ui <- fluidPage(
  titlePanel("Interactive Synteny Analysis Application"),
  sidebarLayout(
    sidebarPanel(
      textInput("gene", "Enter Gene Name:", value = "modA"),
      textInput("product", "Enter Product Name:", value = "alkaline phosphatase D family protein"),
      sliderInput("sidestep", "Select Sidestep (bp):", min = 10000, max = 10000000, value = 1000000),
      sliderInput("window_size", "Select Window Size (bp):", min = 10000, max = 10000000, value = 50000),
      actionButton("plotButton", "Generate Plot")
    ),
    mainPanel(
      plotlyOutput("genePlot")
    )
  )
)

# Server
server <- function(input, output) {
  observeEvent(input$plotButton, {
    # target gene
    target_gene_name <- input$gene
    target_product <- input$product
    sidestep <- input$sidestep
    window_size <- input$window_size
    gtf_folder <- "./data/GTF_files/"
    
    # Получить список файлов GTF в папке
    gtf_files <- list.files(path = gtf_folder, pattern = "\\.gtf$", full.names = TRUE)
    
    # Создать пустой список для хранения результатов
    result_list <- list()
    
    # Цикл по каждому файлу GTF
    for (gtf_file in gtf_files) {
      # Извлечь название штамма из имени файла
      strain_name <- tools::file_path_sans_ext(basename(gtf_file))
      
      # Загрузка GTF файла
      gtf_data <- rtracklayer::readGFF(gtf_file)
      
      # Находим гены в окрестности target_gene_name
      target_region <- gtf_data %>%
        filter(type == "CDS" & gene == target_gene_name) %>%
        mutate(start_target = start, end_target = end)
      
      # Определяем окрестность sidestep оснований в каждую сторону от target
      flanking_region <- target_region %>%
        mutate(start_flanking = start_target - sidestep, end_flanking = end_target + sidestep)
      
      # Фильтруем данные для генов в определенной окрестности
      genes_in_flanking_region <- gtf_data %>%
        filter(type == "CDS" & 
                 end >= flanking_region$start_flanking & start <= flanking_region$end_flanking)
      
      # Проверяем, не является ли набор данных пустым
      if (nrow(genes_in_flanking_region) > 0) {
        # Отбираем гены с продуктом target_product
        target_genes <- genes_in_flanking_region %>%
          filter(product == target_product)
        
        # Добавляем результат к списку, включая название файла
        result_list[[strain_name]] <- target_genes %>%
          mutate(filename = gtf_file)
      } else {
        # Добавляем сообщение к списку, включая название файла
        result_list[[strain_name]] <- data.frame(message = "Genes in the vicinity of target_gene_name not found.", filename = gtf_file)
      }
    }
    
    # Соединяем все результаты в одну таблицу (если необходимо)
    final_result <- bind_rows(result_list, .id = "strain")
    
    # Получение списка файлов GTF в папке
    gtf_files <- final_result$filename
    
    # Функция для центрирования координат относительно указанного гена
    center_coordinates <- function(gtf_data, center_gene, window_size) {
      center_gene_coords <- gtf_data %>%
        filter(type == "CDS", gene_id == center_gene) %>%
        dplyr::select(seqid, center_start = start, center_end = end)
      
      centered_data <- gtf_data %>%
        filter(type == "CDS",
               seqid %in% center_gene_coords$seqid) %>%
        mutate(centered_start = start - center_gene_coords$center_start,
               centered_end = end - center_gene_coords$center_start,
               direction = if_else(centered_start < centered_end, "->", "<-")) %>%
        filter(centered_end >= -window_size, centered_start <= window_size) %>%
        dplyr::select(gene, product, centered_start, centered_end, direction)
      
      return(centered_data)
    }
    
    # Считывание и центрирование данных
    center_gene_list <- final_result$gene_id
    
    centered_genes_list <- lapply(seq_along(gtf_files), function(i) {
      file <- gtf_files[i]
      center_gene <- center_gene_list[i]
      
      gtf <- rtracklayer::readGFF(file)
      centered_genes <- center_coordinates(gtf, center_gene, window_size)
      return(centered_genes)
    })
    
    all_centered_genes <- bind_rows(centered_genes_list, .id = "file_id")
    
    p <- ggplot() +
      geom_rect(data = all_centered_genes,
                aes(xmin = centered_start, xmax = centered_end, ymin = as.numeric(file_id) - 0.4, ymax = as.numeric(file_id) + 0.4, fill = gene),
                alpha = 1) +
      scale_y_continuous(breaks = seq_along(gtf_files), labels = basename(gtf_files), guide = FALSE) +
      theme_minimal() +
      geom_vline(xintercept = 0, color = "blue", linetype = "solid") +
      ggtitle(paste("Target Gene:", target_gene_name, "\n", "Target Product:", target_product)) +
      theme(axis.text.x = element_blank(),  # Удаление подписей оси x
            axis.text.y = element_text(size = 8))  # Уменьшение размера шрифта по вертикали
    
    
    # Преобразование в интерактивный график
    p <- ggplotly(p)
    
    # Вывод графика
    output$genePlot <- renderPlotly({
      return(p)
    })
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
