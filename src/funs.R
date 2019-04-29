
# ---- Setup ----

# Function to load packages
install_packages <- function(install_list, name_list = NULL, 
                             install_cmd = "utils::install.packages", ...) {
  
  # If using devtools::install_github(), the package name will be different 
  # from the repository name in install_list
  if (is.null(name_list)) {
    name_list <- install_list
  }
  
  # Install packages
  for (i in 1:length(install_list)) {
    # require() returns TRUE invisibly if it was able to load package
    if (!require(name_list[i], character.only = T)) {
      
      install_pkg <- strsplit(install_cmd, "::")[[1]][1]
      
      if (!require(install_pkg, character.only = T)) {
        install.packages(install_pkg, dependencies = T, ...)
        require(install_pkg)
      }
      
      install_cmd <- paste0(install_cmd, "(install_list[i], ...)")
      eval(parse(text = install_cmd))
      
      require(name_list[i], character.only = T)
    }
  }
}


# Package lists
Bioc_packages <- c("ComplexHeatmap", "DESeq2")

R_packages <- c(
  "depmixS4",   "ggseqlogo",   
  "gridExtra",  "ggpubr",
  "cowplot",    "rlang",    
  "Rcpp",       "valr",     
  "knitr",      "magrittr", 
  "broom",      "tidyverse", 
  "scales",     "RMySQL",
  "foreach",    "lazyeval",
  "pryr",       "devtools"
)

# Install packages
install_packages(R_packages)

install_packages(
  install_list = Bioc_packages, 
  install_cmd = "BiocManager::install"
)

install_packages(
  install_list = "hadley/multidplyr", 
  name_list    = "multidplyr", 
  install_cmd  = "devtools::install_github"
)

# Directories
proj_dir <- "~/Projects/TFIISdn_paper_old/"
results_dir   <- str_c(proj_dir, "results/")
data_dir      <- str_c(proj_dir, "data/")
gene_list_dir <- str_c(data_dir, "Gene_lists/")
R_list_dir    <- str_c(results_dir,"R_gene_lists/")

# Theme info
theme_info <- theme_classic() +
  theme(
    strip.background  = element_blank(),
    strip.text        = element_text(size = 16, face = "bold"),
    axis.title        = element_text(size = 16, face = "bold"),
    axis.text         = element_text(size = 12, face = "bold", color = "black"),
    axis.line         = element_line(size = 2, lineend = "square"),
    axis.ticks        = element_line(size = 2),
    axis.ticks.length = unit(10, units = "point"),
    legend.title      = element_blank(),
    legend.text       = element_text(size = 12, face = "bold"),
    legend.position   = c(0.8, 0.9),
    legend.background = element_blank()
  )

# Plot colors 
color_list <- list(
  purple_1 = "#8c6bb1",
  purple_2 = "#253494",
  purple_3 = "#80379c",
  orange_1 = "#fe9929",
  orange_2 = "#ec7014",
  orange_3 = "#ff8422",
  brown_1  = "#cc4c02",
  brown_2  = "#993404",
  green_1  = "#78c679",
  green_2  = "#41ab5d",
  green_3  = "#238443",
  green_4  = "#006837",
  green_5  = "#004529",
  blue_1   = "#3288bd",
  blue_2   = "#225ea8",
  blue_3   = "#08519c",
  blue_4   = "#327CD6",
  aqua_1   = "#41b6c4",
  grey_1   = "#969696",
  grey_2   = "#737373",
  grey_3   = "#525252",
  grey_4   = "#252525",
  black    = "black",
  red_1    = "#fc4e2a",
  red_2    = "#e31a1c",
  red_3    = "#bd0026"
)



# ---- General ----

# Function to retrieve colors from color_list
get_color <- function(targets) {
  as.character(color_list[targets])
}

# Function to generate file path for TFIISdn data
get_path <- function(input_file, outer_dir = NULL, 
                     res_dir = results_dir) {

  # Generate name for inner directory containing target file
  inner_dir <- str_split(input_file, "\\.")[[1]][1]
  
  # Generate name for outer data directory
  if (is.null(outer_dir)) {
    outer_dir <- str_split(inner_dir, "_")[[1]][1]
  }
  
  outer_dir <- str_c(outer_dir, "-seq/")
  
  # Output full path
  res <- str_c(res_dir, outer_dir, inner_dir, "/", input_file)
  
  res
}

# Function to read TFIISdn files into nested data.frame
import_dfs <- function(input_list, col_names, get_path = T) {
  
  # input_list must be a list of lists 
  # Each inner list should be comma separated and contain the key name and file name
  # If get_path = T, the data directory name to pass to get_path() can be included as the third value
  
  # Initialize empty data.frame
  res <- tibble(key = character(), data = list())
  
  # Create data.frame containing the key and file path
  for (i in seq_along(input_list)) {
    key_name  <- input_list[[i]][1]
    file_path <- input_list[[i]][2]
    
    # Use get_path()
    if (get_path) {
      if (length(input_list[[i]]) < 3) {
        file_path <- get_path(file_path)
        
      } else {
        dir_name  <- input_list[[i]][3]
        file_path <- get_path(file_path, dir_name)
      }
    }
    
    # Create data.frame for new line
    new_row <- tibble(key = key_name, data = list(file_path))
    
    # Add new line to output data.frame
    res <- bind_rows(res, new_row)
  }
  
  # Import files into data.frame
  res <- res %>%
    mutate(data = map(data, read_tsv, col_names))
  
  res
}

# Function to add key column to each file in list
add_key <- function(input, rename_col = NULL) {
  
  # If input file is a data.frame
  if (is_tibble(input)) {
    
    # If a column should be renamed
    if (!is.null(rename_col)) {
      old_col <- sym(rename_col)
      
      res <- map2(input$data, input$key, ~rename(.x, !!sym(.y) := !!old_col)) %>%  # Had trouble getting this to work using dplyr::mutate()
        tibble(key = input$key, data = .)
      
    # If a new column should be added  
    } else {
      res <- input %>%
        mutate(data = map2(data, key, ~mutate(.x, key = .y)))
    }
    
  # If the input file is a named list 
  } else {
    
    # If a column should be renamed
    if (!is.null(rename_col)) {
      res <- imap(input, ~{
        .x %>%
          set_names(str_replace(colnames(.x), rename_col, .y))
      })
      
    # If a new column should be added
    } else {
      res <- imap(input, ~mutate(.x, key = .y))
    }
  }
  
  res
  
  
  # OLD VERSION
  # add_key <- function(input, col_name = "count", gather_keys = F) {
  #   
  #   col_sym <- sym(col_name)
  #   
  #   # If input file is a data.frame
  #   if (is_tibble(input)) {
  #     res <- map2(input$data, input$key, ~ {
  #       new_name <- sym(.y)
  #       
  #       out_df <- .x %>% 
  #         dplyr::rename(!!new_name := !!col_sym)
  #       
  #       if (gather_keys) {
  #         out_df <- out_df %>%
  #           gather(key, !!col_sym, !!new_name)
  #       }
  #       
  #       out_df
  #     })
  #     
  #   # If input file is a named list
  #   } else {
  #     res <- imap(input, ~ {
  #       out_df <- .x
  #       
  #       colnames(out_df) <- str_replace(colnames(.x), col_name, .y)
  #       
  #       if (gather_keys) {
  #         out_df <- out_df %>%
  #           gather(key, !!col_sym, .y)
  #       }
  #       
  #       out_df
  #     })
  #   }
  #   
  #   res
  #   
  # dir(gene_list_dir, "hg19\\w*ksep.bed|hg19_1klong.bed")
  # dir("\\.csv$") %>%
  #   set_names() %>%
  #   map(read.csv) %>%
  #   imap(~ transform(.x, filename = .y))
  # }
}

# Function to generate key names for plots
rename_key <- function(input_df) {
  
  # Function to find and replace strings in key
  search_replace <- function(in_str, in_search, in_replace) {
    
    search_str <- str_c(".*(", in_search, ")\\w*")
    
    res <- ifelse(
      grepl(in_search, in_str), 
      str_replace(in_str, search_str, in_replace), 
      in_str
    )
    
    res
  }
  
  res <- input_df
  
  # Only add sense and anti-sense labels if both are present in data.frame
  if (any(grepl("_sen$", input_df$key)) && any(grepl("_anti$", input_df$key))) {
    res <- res %>%
      mutate(
        key = ifelse(grepl("_anti$", key), str_c(key, " anti-sense"), key),
        key = ifelse(grepl("_sen$", key),  str_c(key, " sense"), key)
      )
  }
  
  # Rename keys
  res <- res %>%
    mutate(
      key = search_replace(key, "H2O", "Control"),
      key = search_replace(key, "Dox|IISdn", "+TFIISDN"),
      key = search_replace(key, "IISwt", "+TFIISWT")
    )
  
  res
}

# Function to retrieve unique gene names
get_uniq_names <- function(input_df) {
  input_df %>%
    dplyr::select(name) %>%
    unique()
}

# Function to merge data.frames 
merge_dfs <- function(input_df, na_omit = T, merge_by = "name") {
  
  res <- input_df %>%
    add_key(rename_col = "count")
  
  res <- res$data %>%
    reduce(full_join, by = merge_by)
  
  if (na_omit) {
    res <- na.omit(res)
  }
  
  res <- res %>%
    gather(key, count, -name)
  
  res
}

# Function to merge metaplot bed files 
merge_wins <- function(input_df, win_num, ref_win = 1, dist_col = "win_id",
                       win_min = 1, win_max = win_num, ignore_anti = F) {
  
  # Merge bed files, print memory usage and system time
  system.time({
    mem_change({
      
      dist_sym <- sym(dist_col)
      res      <- ungroup(input_df)
      
      # Filter nested bed file data.frame
      res <- res %>%
        mutate(data = map2(data, key, ~{
          .x <- .x %>%
            group_by(name) %>%
            filter(
              # min(win_id) == 1,
              # max(win_id) == eval(win_num),
              n()    == eval(win_num),
              win_id >= eval(win_min),
              win_id <= eval(win_max)
            )
          
          # Remove genes with zero counts
          # Ignore anti-sense data when filtering based on signal
          if (!ignore_anti | !grepl("anti", .y)) {
            .x <- .x %>%
              filter(sum(count) > 0)
          }

          .x
        }))
      
      # Create list of overlapping genes
      shared_genes <- res$data %>%
        map(get_uniq_names) %>%
        purrr::reduce(semi_join, by = "name")
      
      # Merge tables
      res <- res %>%
        ungroup() %>%
        mutate(data = map2(data, key, ~{
          .x <- .x %>%
            ungroup() %>%
            semi_join(shared_genes, by = "name")

          # Express win_id as kb from reference window
          if (!is.null(ref_win)) {
            ref_win <- eval(ref_win)
            len     <- unique(.x$end - .x$start)

            if (length(len) == 1) {
              .x <- .x %>%
                mutate(
                  !!dist_sym := (win_id - ref_win) * len / 1000,
                  win_len     = len
                ) %>%
                arrange(name, !!dist_sym)

            } else {
              .x <- .x %>%
                mutate(win_len = (end - start) / 1000) %>%
                arrange(name, win_id) %>%
                group_by(name) %>%
                mutate(
                  id          = win_id,
                  !!dist_sym := cumsum(win_len),
                  !!dist_sym := lag(!!dist_sym, 1, default = 0),
                  ref_dist    = ifelse(id == ref_win, !!dist_sym, 0),
                  ref_dist    = max(ref_dist),
                  !!dist_sym := !!dist_sym - ref_dist
                ) %>%
                ungroup() %>%
                dplyr::select(-ref_dist, -id)
            }
          }

          .x <- .x %>%
            dplyr::select(-chrom, -start, -end, -strand)

          .x
        })) %>%
        add_key()
      
    }) %>%
      cat(sep = "", "Memory usage:\n  ", ., " B\n\nSystem time:\n")
    
  }) %>%
    print()
  
  res
}

# Function to calculate mean signal
calc_mean_signal <- function(input_df, gene_list = NULL, data_col = "count", 
                             group_cols = c("type", "key", "rep", "strand"), 
                             split_key = T, rel_freq = F, rep_mean = F) { 
  
  # Calculate mean signal, print memory usage and system time
  system.time({
    mem_change({
      
      res <- input_df %>% 
        ungroup()
      
      # Grouping columns
      data_sym   <- sym(data_col)
      group_syms <- syms(group_cols)
      
      # Unnest nested data.frame
      if (all(c("key", "data") %in% names(res))) {
        res <- res %>%
          dplyr::select(data) %>%
          unnest()
      }
      
      # Merge data.frame() with gene_list
      if (!is.null(gene_list)) {
        res <- res %>%
          semi_join(gene_list, by = "name")
      }
      
      # Split key column into grouping columns
      if (split_key) {
        res <- res %>%
          separate(key, sep = "_", into = group_cols)
      }
      
      # Calculate relative signal ----
      if (rel_freq) {
        # By default strand is ignored when calculating relative signal
        no_strand_cols <- grep("strand", group_cols, value = T, invert = T)
        no_strand_syms <- syms(no_strand_cols)
        
        res <- res %>%
          group_by(!!!no_strand_syms, name) %>%
          mutate(!!data_sym := !!data_sym / sum(!!data_sym)) %>%
          ungroup()
      }
      
      # Calculate mean signal ----
      res <- res %>%
        group_by(!!!group_syms) %>%
        mutate(n_gene = n_distinct(name)) %>% 
        group_by(win_id, n_gene, add = T) %>% 
        summarize(!!data_sym := mean(!!data_sym)) %>%
        ungroup()
      
      # Calculate mean and SEM for replicates ----
      if (rep_mean) {
        no_rep_cols <- grep("rep", group_cols, value = T, invert = T)
        group_syms <- syms(no_rep_cols)
        
        res <- res %>%
          group_by(!!!group_syms) %>%
          mutate(n_rep  = n_distinct(rep)) %>%
          group_by(win_id, n_rep, n_gene, add = T) %>%
          summarize(
            SEM_count   = sd(!!data_sym),
            !!data_sym := mean(!!data_sym)
          ) %>%
          ungroup() %>%
          mutate(SEM_count = SEM_count / sqrt(n_rep))
      }
      
      # Reunite grouping columns
      if (split_key) {
        res <- res %>%
          unite(key, !!!group_syms, sep = "_")
      }
      
    }) %>%
      cat(sep = "", "Memory usage:\n  ", ., " B\n\nSystem time:\n")
    
  }) %>%
    print()
  
  res
}

# Function to add comma to gene number
add_commas <- function(num) {
  
  nums <- str_split(num, "")[[1]]
  nums_len <- length(nums)
  nums_rev <- rev(nums)
  
  res <- imap(nums_rev, ~ {
    if (.y %% 3 == 0 & .y != nums_len) {
      .x <- c(.x, ",")
    }
    
    .x
  }) %>%
    unlist() %>%
    rev() %>%
    str_c(collapse = "")
  
  res
}

# Function to add subscripts to plot labels
add_subscript <- function(input_str, pattern = "DN|WT", bold = T) {
  
  # Extract characters and position of pattern
  pos     <- str_locate(input_str, pattern)
  str_len <- str_length(input_str)
  chars   <- str_split(input_str, "") %>% 
    unlist()
  
  # Convert pattern to subscript
  if (str_detect(input_str, pattern)) {
    start_chars <- ""
    pattern     <- chars[pos[1] : pos[2]] %>% 
      str_c(collapse = "")
    
    if (pos[1] > 1) {
      start_chars <- chars[1 : (pos[1] - 1)] %>%
        str_c(collapse = "")
    }
    
    end_chars <- ""
    
    if (pos[2] < str_len) {
      end_chars <- chars[(pos[2] + 1) : str_len] %>%
        str_c(collapse = "")
    }
    
    res <- bquote(.(start_chars) * ""[.(pattern)] * .(end_chars))
    
    if (bold) {
      res <- bquote(bold(.(start_chars) * ""[.(pattern)] * .(end_chars)))
    }
    
  } else {
    start_chars <- chars[1 : str_len] %>%  # This was only way to get facet_wrap labeller working
      str_c(collapse = "")
    
    res <- start_chars
    
    if (bold) {
      res <- bquote(bold(.(start_chars)))
    }
  }
  
  res
}

# Function to create metaplots 
create_metaplots <- function(input_df, data_col = "count", stranded = F, 
                             plot_SEM = F, smooth_line = F, trans = 1,
                             plot_levels = NULL, plot_colors = NULL, 
                             x_lim = NULL, x_breaks = NULL, x_labs = NULL, 
                             y_lim = NULL, y_breaks = NULL, n_label = c(0.81, 0.83),
                             v_line = NULL, h_line = NULL, ...) {
  
  # Modify input data.frame ----
  
  # Identify data column
  data_sym  <- sym(data_col)
  
  # Transform counts
  sen_data <- input_df %>% 
    ungroup() %>%
    mutate(
      !!data_sym := !!data_sym * trans,
      key         = fct_relevel(key, plot_levels)
    )
  
  # Transform SEM
  if (plot_SEM) {
    SEM_col <- names(input_df) %>%
      grep("SEM", ., value = T) %>%
      sym()
    
    sen_data <- sen_data %>%
      mutate(!!SEM_col := !!SEM_col * trans)
  }
  
  # Adjust x-axis limits
  if (!is.null(x_lim)) {
    sen_data <- sen_data %>%
      filter(
        win_id >= x_lim[1],
        win_id <= x_lim[2]
      )
  }
  
  # Identify anti-sense data
  if (stranded) {
    
    # Function to split data based on strand
    get_strand_data <- function(input_df, target_str, flip = F) {
      res <- input_df %>%
        filter(grepl(target_str, key)) %>%
        mutate(
          key = str_remove(key, target_str),
          key = fct_relevel(key, plot_levels)
        )
      
      if (flip) {
        res <- res %>% 
          mutate(!!data_sym := !!data_sym * -1)
      }
      
      res
    }
    
    plot_levels <- str_remove(plot_levels, " anti-sense| sense") %>%
      unique()
    
    anti_data <- get_strand_data(sen_data, " anti-sense", flip = T)
    sen_data  <- get_strand_data(sen_data, " sense")
  }
  
  
  # Create ribbon for standard error ----
  
  # Create ggplot object
  meta_plot <- sen_data %>% 
    ggplot(aes(win_id, !!data_sym, color = key, fill = key))
  
  # Plot ribbon for SEM 
  if (plot_SEM) {
    
    # Function to plot ribbon
    plot_ribbon <- function(input_df) {
      geom_ribbon(
        aes(
          ymin     = !!data_sym - !!SEM_col,
          ymax     = !!data_sym + !!SEM_col,
          linetype = NA,
          fill     = key
        ),
        data        = input_df, 
        show.legend = F,
        alpha       = 0.25
      )
    }
    
    meta_plot <- meta_plot +
      plot_ribbon(sen_data)
    
    if (stranded) {
      meta_plot <- meta_plot +
        plot_ribbon(anti_data)
    }
  }
  
  
  # Plot lines for mean signal ----
  
  # Plot non-faded lines for sense signal
  if (smooth_line) {
    meta_plot <- meta_plot +
      geom_line(
        data   = sen_data, 
        stat   = "smooth", 
        method = "loess", 
        se     = F, 
        ...
      )
    
  } else {
    meta_plot <- meta_plot +
      geom_line(data = sen_data, ...)      
  }
  
  # Plot faded lines for antisense signal
  if (stranded) {
    
    # Function to generate antisense lines
    plot_anti_line <- function(...) {
      
      if (smooth_line) {
        res <- geom_line(
          data   = anti_layer, 
          stat   = "smooth", 
          method = "loess", 
          se     = F, 
          ...
        )
        
      } else {
        res <- geom_line(data = anti_layer, ...)
      }
      
      res
    }
    
    # Plot faded lines for mean signal
    for (i in seq_along(plot_levels)) {
      anti_layer <- anti_data %>%
        filter(key == plot_levels[i])
      
      plot_color <- plot_colors[i]
      
      meta_plot <- meta_plot +
        plot_anti_line(color = "white", ...) +
        plot_anti_line(
          color = plot_color, 
          alpha = 0.5,
          ...
        )
    }
  }
  
  
  # Add N label ----
  
  if (!is.null(n_label)) {
    x_min <- min(sen_data$win_id)
    x_max <- max(sen_data$win_id)
    y_min <- min(sen_data[, data_col])
    y_max <- max(sen_data[, data_col])
    
    if (stranded) {
      y_min <- min(anti_data[, data_col])
    }
    
    # Account for SEM when calculating y-axis limits
    if (plot_SEM) {
      plot_limits <- sen_data %>%
        mutate(
          min_count = count - !!SEM_col,
          max_count = count + !!SEM_col
        )
      
      y_min <- min(plot_limits$min_count)
      y_max <- max(plot_limits$max_count)
      
      if (stranded) {
        plot_limits <- anti_data %>%
          mutate(min_count = count - !!SEM_col)
        
        y_min <- min(plot_limits$min_count)
      }
    }
    
    # Adjust label position
    if (!is.null(x_lim)) {
      x_min <- x_lim[1]
      x_max <- x_lim[2]
    }
    
    if (!is.null(y_lim)) {
      y_min <- y_lim[1]
      y_max <- y_lim[2]
    }
    
    n     <- unique(input_df$n_gene)
    n_lab <- str_c("N = ", add_commas(n))
    n_x   <- x_min + ((x_max - x_min) * n_label[1])
    n_y   <- y_min + ((y_max - y_min) * n_label[2])
    
    # Add N label to plot
    meta_plot <- meta_plot +
      geom_text(
        aes(x = n_x, y = n_y, label = n_lab),
        color         = "black",
        fontface      = "bold",
        check_overlap = T
      )
  }
  
  
  # Options for plot aesthetics ----
  
  # Set plot colors
  if (!is.null(plot_colors)) {
    plot_labels <- map(plot_levels, add_subscript)
    
    meta_plot <- meta_plot +
      scale_color_manual(
        values = plot_colors, 
        labels = plot_labels
      ) +
      scale_fill_manual(
        values = plot_colors, 
        labels = plot_labels
      )
  }
  
  # Adjust x-breaks
  if (!is.null(x_breaks)) {
    meta_plot <- meta_plot +
      scale_x_continuous(
        breaks = x_breaks, 
        labels = x_labs
      )
  }
  
  # Adjust y-breaks
  if (!is.null(y_breaks)) {
    meta_plot <- meta_plot +
      scale_y_continuous(breaks = y_breaks)
  }
  
  # Add vertical lines
  if (!is.null(v_line)) {
    meta_plot <- meta_plot + 
      geom_vline(
        xintercept = v_line, 
        size       = 1, 
        linetype   = 2
      )
  }
  
  # Add horizontal lines
  if (!is.null(h_line)) {
    meta_plot <- meta_plot +
      geom_hline(
        yintercept = h_line, 
        size       = 1
      )
  }
  
  # Add theme info
  meta_plot <- meta_plot + 
    coord_cartesian(ylim = y_lim) +
    theme_info +
    theme(axis.title.x = element_blank())
  
  meta_plot
}

# Function to plot elongation rates 
create_boxplots <- function(input_df, data_col = "count", log_trans = F, 
                            plot_levels = NULL, plot_colors = NULL, ticks = T,
                            y_lim = NULL, y_breaks = NULL, n_groups = "key", 
                            n_label = c(2, 0.001), h_line = NULL, ...) {
  
  # Create boxplots ----
  
  # Identify data column
  data_sym <- sym(data_col)
  
  # Set plot levels
  input_df <- input_df %>%
    ungroup() %>%
    mutate(key = fct_relevel(key, plot_levels))
  
  # Log transform data
  if (log_trans == T) {
    input_df <- input_df %>%
      mutate(!!data_sym := log2(!!data_sym))
  }
  
  # Add columns for N label
  if (!is.null(n_label)) {
    
    if (!is.null(n_groups)) {
      n_groups <- syms(n_groups)
      
      input_df <- input_df %>%
        group_by(!!!n_groups)
    }
    
    input_df <- input_df %>%
      mutate(
        n_gene = n_distinct(name),
        y_min  = min(!!data_sym),
        y_max  = max(!!data_sym)
      ) %>%
      ungroup() %>%
      nest(n_gene, .key = n_gene) %>%
      mutate(n_lab = map(n_gene, ~str_c("N = ", add_commas(.x)))) %>%  # Not sure how to make this work without nest()
      unnest()
    
    if (!is.null(y_lim)) {
      input_df <- input_df %>%
        mutate(
          y_min = y_lim[1],
          y_max = y_lim[2]
        )
    }
  }
  
  # Create boxplot
  box_plot <- input_df %>% 
    ggplot(aes(key, !!data_sym, color = key)) +
    geom_boxplot(...) +
    coord_cartesian(ylim = y_lim) +
    theme_info +
    theme(
      legend.position = "none",
      axis.title.x = element_blank()
    ) +
    scale_x_discrete(labels = map(plot_levels, add_subscript))
  
  
  
  # Options for plot aesthetics ---- 
  
  # Set plot colors
  if (!is.null(plot_colors)) {
    plot_labels <- map(plot_levels, add_subscript)
    
    box_plot <- box_plot +
      scale_color_manual(
        values = plot_colors, 
        labels = plot_labels
      )
  }
  
  # Remove tick marks
  if (ticks == F) {
    box_plot <- box_plot + 
      theme(axis.ticks.x = element_blank())
  }
  
  # Adjust y-breaks
  if (!is.null(y_breaks)) {
    box_plot <- box_plot + 
      scale_y_continuous(breaks = y_breaks)
  }
  
  # Add horizontal lines
  if (!is.null(h_line)) {
    box_plot <- box_plot +
      geom_hline(
        yintercept = h_line, 
        size       = 0.75, 
        linetype   = 1
      )
  }
  
  # Add N label
  if (!is.null(n_label)) {
    box_plot <- box_plot +
      geom_text(
        aes(
          x     = n_label[1], 
          y     = y_min + ((y_max - y_min) * n_label[2]),
          label = n_lab
        ),
        color         = "black",
        fontface      = "bold",
        hjust         = 0,
        check_overlap = T
      )
  }
  
  box_plot
}

# Function to generate genome browser tracks
create_browser_tracks <- function(bg_list, gene_name, omit_gene = ",MIR|,SNOR",
                                  both = 0, left = 0, right = 0, win_size,
                                  step_size = NULL, rev_track = F, overlay = F,
                                  y_min = NULL, y_max = NULL, y_scale = "strand_equal",
                                  v_line = "TSS", v_colors = "black", rev_legend = F,
                                  highlight_region = NULL, highlight_color = NULL,
                                  theme_opts = NULL, gene_color = "black",
                                  gene_label = T, arrow_color = "#e31a1c",
                                  overlap_genes = F, track_heights = NULL,
                                  dist_scale = c(0.7, 0.9), ...) {
  
  # Theme info
  browser_theme <- theme_classic() +
    theme_info +
    theme(
      strip.text   = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 14, color = "black", face = "bold"),
      axis.text.x  = element_blank(),
      axis.text.y  = element_text(size = 10, color = "black", face = "bold")
    )
  
  
  # Create windows ----
  
  # Identify target region
  gene_coords <- gene_annotations %>%
    dplyr::select(-data) %>% 
    filter(grepl(gene_name, name))
  
  if (nrow(gene_coords) > 1) {
    return(gene_coords$name)
  }
  
  # Identify the TSS and pAS
  strand <- gene_coords$strand
  
  if (strand == "+") {
    TSS_coord <- gene_coords$start
    pAS_coord <- gene_coords$end
    
  } else {
    TSS_coord <- gene_coords$end
    pAS_coord <- gene_coords$start
  }
  
  # Adjust region coordinates
  if (both > 0) {
    left  <- both
    right <- both
  }
  
  # Identify target region
  target_region <- gene_coords %>%
    mutate(
      start = ifelse(
        rev_track, 
        start - right, 
        start - left
      ),
      end = ifelse(
        rev_track, 
        end + left, 
        end + right
      )
    )
  
  region_start <- target_region$start
  region_end   <- target_region$end
  
  # Create windows
  if (!is.null(step_size)) {
    
    if (step_size < win_size) {
      target_wins <- target_region %>%
        bed_makewindows(
          win_size  = win_size, 
          step_size = step_size
        )
      
    } else {
      warning("Step size must be smaller than window size")
      
      target_wins <- target_region %>% 
        bed_makewindows(win_size = win_size)
    }
    
  } else {
    target_wins <- target_region %>% 
      bed_makewindows(win_size = win_size)
  }
  
  
  # Merge bedGraphs ----
  
  bg_merge <- imap(bg_list, ~ {
    
    # Extract track options
    opt_cols <- .y %>%
      str_split(",") %>%
      map(str_trim) %>%
      unlist() %>%
      str_trim()
    
    # Add new columns for bedGraph options
    flip <- F
    
    new_bg <- .x %>%
      mutate(
        group      = opt_cols[1],
        group_id   = as.numeric(opt_cols[2]),
        treatement = opt_cols[3],
        color      = opt_cols[4],
        flip_track = flip
      )
    
    # Plot signal as negative numbers
    if ("flip" %in% opt_cols) {
      flip <- T
      
      new_bg <- new_bg %>%
        mutate(
          count      = count * -1,
          flip_track = flip
        )
    }
    
    # data.frame for regions containing signal
    signal <- target_wins %>%
      bed_intersect(y = new_bg) %>%
      mutate(count = count.y * .overlap) %>% 
      dplyr::select(
        chrom,           start    = 2,  
        end        = 3,  .win_id  = 7, 
        group      = 11, group_id = 12,
        treatment  = 13, color    = 14,
        flip_track = 15, count
      ) %>% 
      group_by(.win_id) %>%
      mutate(count = sum(count)) %>% 
      ungroup() %>% 
      unique() %>%
      mutate(count = count / (end - start))   
    
    # data.frame for regions lacking signal
    no_signal <- target_wins %>% 
      dplyr::select(chrom, start, end, .win_id) %>%
      bed_subtract(signal) %>%
      mutate(
        group      = opt_cols[1],
        group_id   = as.numeric(opt_cols[2]),
        treatment  = opt_cols[3],
        color      = opt_cols[4],
        flip_track = flip,
        count      = 0
      ) 
    
    bind_rows(signal, no_signal)
  }) %>%
    bind_rows()
  
  
  # Adjust coords and identify min and max values ----
  
  # Function to adjust bedGraph coordinates
  adjust_coords <- function(input_df) {
    
    if (is.null(step_size)) {
      step_size <- -1
    }
    
    cols_to_keep <- c(
      "group",      "treatment", 
      "group_id",   "color", 
      "flip_track", "count"
    )
    
    start_coords <- input_df %>%
      dplyr::select(coord = 2, cols_to_keep)
    
    end_coords <- input_df %>%
      mutate(end = start + win_size - step_size - 1) %>% 
      dplyr::select(coord = 3, cols_to_keep)
    
    # Identify min and max y values
    res <- bind_rows(start_coords, end_coords) %>%
      mutate(
        plot_min = min(count),
        plot_max = max(count)
      ) %>%
      group_by(group_id) %>%
      mutate(
        group_min = min(count),
        group_max = max(count)
      ) %>%
      ungroup()
    
    res
  }
  
  signal_df <- adjust_coords(bg_merge)
  
  # Create data.frame for highlight_region
  if (!is.null(highlight_region)) {
    
    # if (is.tibble(highlight_region)) {
    #   highlight_region <- list(highlight_region)
    # }
    
    highlight_dfs <- highlight_region %>%
      map(~{
        .x %>%
          bed_intersect(x = bg_merge, suffix = c("", ".y")) %>%
          filter(
            .overlap > 0,
            count    > 0
          ) %>%
          adjust_coords()
      })
  }
  
  # Set position of distance scale
  max_group   <- max(signal_df$group_id)
  x_min       <- min(signal_df$coord)
  x_max       <- max(signal_df$coord)
  x_range     <- x_max - x_min
  
  if (!is.null(dist_scale)) {
    scale_range <- round(x_range * 0.05)
    scale_start <- x_min + (x_range * dist_scale[1])
    scale_end   <- scale_start + scale_range
    scale_lab   <- round(scale_range / 1000, digits = 1) %>%
      str_c("  ", .," kb  ")
    
    if (rev_track) {
      scale_start <- x_max - (x_range * dist_scale[1])
      scale_end   <- scale_start - scale_range
    }
  }
  
  
  # Create signal tracks ---- 
  
  # Function to create signal track
  create_signal_track <- function(input_df, input_id, ...) {
    
    # Format y-axis scale ----
    if (is.null(y_min)) {
      y_min <- unique(input_df$group_min) %>%
        round(digits = 1)
        # This ifelse statement causes some tracks to appear blank when using
        # certain window sizes. This bug was introduced in commit a6e1da2 
        # on Dec 12, 2018.
        # round(digits = ifelse(. >= 0.1, 1, 2))
    }
    
    if (is.null(y_max)) {
      y_max <- unique(input_df$group_max) %>%
        round(digits = 1)
        # This ifelse statement causes some tracks to appear blank when using
        # certain window sizes. This bug was introduced in commit a6e1da2 
        # on Dec 12, 2018.
        # round(digits = ifelse(. >= 0.1, 1, 2))
    }
    
    # Set negative and positive y limits equal
    if (y_scale == "strand_equal" & y_min < 0) {
      max_value <- max(abs(y_min), abs(y_max))
      
      y_min <- max_value * -1
      y_max <- max_value
    }
    
    y_lims   <- c(y_min, y_max)
    y_breaks <- c(0, y_max)
    
    if (y_min < 0) {
      y_breaks <- c(y_min, y_breaks)
    }
    
    # Retrieve plot levels, colors, and labels ----
    plot_colors <- input_df %>%
      dplyr::select(treatment, color) %>%
      unique() %>%
      dplyr::select(color) %>%
      unlist() %>%
      as.character()
    
    label_aes <- input_df %>%
      filter(!flip_track) %>%
      dplyr::select(color, treatment) %>%
      unique()
    
    plot_labels <- label_aes$treatment %>%
      map(add_subscript)
    
    label_colors     <- label_aes$color
    label_n          <- length(plot_labels)
    treatment_levels <- unique(input_df$treatment)
    treatment_n      <- length(treatment_levels)
    group_levels     <- unique(input_df$group)
    y_label          <- add_subscript(group_levels[1])
    
    # Create signal tracks ----
    res <- input_df %>%
      mutate(
        group     = fct_relevel(group, group_levels),
        treatment = fct_relevel(treatment, treatment_levels)
      ) %>%
      ggplot(aes(coord, count, fill = treatment, color = treatment)) +
      geom_line(size = 0) +  # Need this layer to include the "line" shape in legend 
      geom_ribbon(
        aes(ymin = 0, ymax = count), 
        show.legend = F, 
        ...
      ) +
      scale_y_continuous(
        breaks = y_breaks,
        limits = y_lims
      )
    
    # Add layer for highlighted region
    if (!is.null(highlight_region)) {
      
      highlight_data <- highlight_dfs[[input_id]] %>%
        filter(group_id == input_id)
      
      res <- res +
        geom_point(
          data        = highlight_data,
          color       = highlight_color,
          size        = 0.3,
          show.legend = F
        )
            
      # for (i in seq_along(highlight_groups)) {
      #   res <- res +
      #     geom_point(
      #       data        = highlight_groups[[i]],
      #       color       = highlight_color[i],
      #       size        = 0.15,
      #       show.legend = F
      #     )
      # }
    }
    
    # Modify legend ----
    
    # Modify size, shape and color of legend symbols
    legend_args <- list(
      size     = c(rep(2, label_n)),
      linetype = c(rep("solid", label_n)),
      color    = label_colors 
    )
    
    # Do not want to show legend entries for antisense signal
    if (any(input_df$flip_track)) {
      # Need these NAs to remove blank legend entries
      legend_args$color <- c(label_colors, rep(NA, label_n))
    }
    
    guide_args <- guide_legend(override.aes = legend_args)
    
    if (rev_legend) {
      legend_args$color <- rev(label_colors)
      guide_args        <- guide_legend(reverse = T, override.aes = legend_args)
    }
    
    res <- res +
      scale_color_manual(
        values = plot_colors, 
        labels = plot_labels,
        guide  = guide_args
      )
    
    # Adjust theme elements ----
    res <- res +
      scale_fill_manual(values = plot_colors) +
      browser_theme +
      theme(
        legend.position      = c(0.2, 0.9),
        legend.justification = c(0, 1),
        panel.spacing        = unit(0.5, "cm"), 
        axis.line.x          = element_blank(),
        axis.ticks.x         = element_blank(),
        plot.margin          = unit(c(0.1, 0.2, 0, 0.2), "cm")
      ) +
      coord_cartesian(
        xlim = c(target_region$start, target_region$end)
      ) +
      labs(y = y_label)
    
    # Overlay plots
    if (!overlay) {
      res <- res + 
        facet_wrap(~treatment, ncol = 1)
      
    } else {
      label_n <- 1  # Overlaid plots should only have one distance scale
    }
    
    # Reverse x-axis
    if (rev_track) {
      res <- res + 
        scale_x_reverse()
    }
    
    # Add vertical lines relative to TSS and pAS positions
    if (!is.null(v_line)) {
      
      # A list of vectors can be given to modify groups separately
      if (typeof(v_line) == "list") {
        v_line <- v_line[[input_id]]
      }
      
      if (typeof(v_colors) == "list") {
        v_colors <- v_colors[[input_id]]
      }
      
      # Convert v_line inputs to genomic coordinates
      other_coords <- v_line %>%
        grep("TSS|pAS", ., value = T, invert = T) %>%
        as.numeric() + TSS_coord
        # as.numeric() + pAS_coord
      
      # Add TSS and/or pAS coordinates
      v_line <- v_line %>%
        grep("TSS|pAS", ., value = T) %>%
        str_replace("TSS", as.character(TSS_coord)) %>%
        str_replace("pAS", as.character(pAS_coord)) %>%
        as.numeric() %>%
        c(other_coords)
      
      res <- res +
        geom_vline(
          xintercept = v_line, 
          color      = v_colors, 
          size       = 0.5, 
          linetype   = 2
        )
    }
    
    # Add distance scale to first track
    if (input_id == 1 & !is.null(dist_scale)) {
      
      # Only make the scale visible on upper panel
      scale_colors <- c("black", rep(NA, label_n - 1))
      
      # Add distance scale
      res <- res +
        annotate(
          geom  = "segment",
          size  = 1,
          x     = scale_start,
          xend  = scale_end,
          y     = y_max * dist_scale[2], 
          yend  = y_max * dist_scale[2],
          color = scale_colors
        ) +
        annotate(
          geom     = "text",
          label    = scale_lab,
          size     = 3.5,
          fontface = 2,
          x        = scale_start,
          y        = y_max * dist_scale[2],
          hjust    = 1,
          color    = scale_colors
        )
    }
    
    # Add bottom line
    if (input_id == max_group) {
      res <- res +
        theme(axis.line.x = element_line(size = 2))
    }
    
    # Additional theme adjustments can be passed as a list
    if (!is.null(theme_opts)) {
      
      # A named list of lists can be given to modify groups separately
      if (!is.null(names(theme_opts))) {
        theme_id   <- str_c("group_", input_id)
        theme_opts <- eval(theme_opts[[theme_id]])
      }
      
      res <- res +
        eval(theme_opts)
    }
    
    res
  }
  
  signal_tracks <- signal_df %>%
    nest(-group_id) %>%
    mutate(data = map2(data, group_id, create_signal_track)) %>%
    arrange(group_id)
  
  
  # Create Refseq gene tracks ----
  
  # Identify genes within target region
  target_genes <- target_region %>%
    bed_intersect(gene_annotations, .) %>%
    filter(!grepl(omit_gene, name.x))
  
  target_genes %>% 
    select(chrom, start = 2, end = 3, name = 4) %>%
    print()
  
  target_genes <- target_genes %>%
    dplyr::select(data.x) %>%
    unnest()
  
  # Create gene annotation track
  refseq_tracks <- target_genes %>%
    ggplot() +
    geom_rect(
      aes(
        xmin = start, 
        xmax = end, 
        ymin = y_min, 
        ymax = y_max
      ),
      color = gene_color, 
      fill  = gene_color
    ) +
    scale_y_continuous(breaks = 4) +
    browser_theme +
    theme(
      legend.position = "none",
      axis.line       = element_blank(),
      axis.ticks      = element_blank(),
      axis.title      = element_blank(),
      axis.text       = element_blank(),
      plot.margin     = unit(c(0, 0.2, 0.2, 0.2), "cm")
    ) +
    coord_cartesian(
      xlim = c(region_start, region_end),
      ylim = c(0, 4.75)
    )
  
  # Overlap gene tracks
  if (!overlap_genes) {
    refseq_tracks <- refseq_tracks +
      facet_wrap(~name, ncol = 1)
  }
  
  # Reverse x-axis
  if (rev_track) {
    refseq_tracks <- refseq_tracks +
      scale_x_reverse()
  }
  
  
  # Add arrow to Refseq gene tracks ----
  
  if (gene_label) {
    
    # Adjust arrow so it is in line with the TSS
    if (TSS_coord > region_start & TSS_coord < region_end) {
      arrow_start <- TSS_coord
      
    } else {
      arrow_start <- region_start
      
      if (rev_track) {
        arrow_start <- region_end
      }
    }
    
    # Adjust arrow so it is in line with the pAS
    if (pAS_coord > region_start & pAS_coord < region_end) {
      arrow_end <- pAS_coord
      
    } else {
      arrow_end <- region_end
      
      if (rev_track) {
        arrow_end <- region_start
      }
    }
    
    # Adjust gene label position
    label_x    <- TSS_coord - (x_range * 0.03)
    label_just <- 1
    
    if (rev_track) {
      label_x    <- TSS_coord + (x_range * 0.03)
      label_just <- 1
    }
    
    # Add gene label and arrow to refseq tracks
    gene_name <- str_remove(gene_name, "\\$")
    
    if (grepl(",", gene_name)) {
      gene_name <- str_split(gene_name, ",")[[1]][2] 
    }
    
    refseq_tracks <- refseq_tracks +
      geom_text(
        aes(label_x, 2, label = gene_name), 
        size          = 4,
        fontface      = "bold",
        hjust         = label_just,
        check_overlap = T
      ) +
      geom_segment(
        aes(x = arrow_start, xend = arrow_end, y = 1.5, yend = 1.5),
        size     = 2,
        color    = arrow_color,
        linejoin = "mitre",
        arrow    = arrow(
          type   = "closed", 
          angle  = 30, 
          length = unit(0.1, "pt")
        )
      )
  }
  
  
  
  # Combine signal and refseq tracks ----
  
  browser_tracks <- signal_tracks %>% 
    add_row(
      group_id = "gene", 
      data     = list(refseq_tracks)
    )
  
  if (is.null(track_heights)) {
    track_heights <- c(rep(1, nrow(signal_tracks)), 0.17)
  }
  
  final_plots <- plot_grid(
    plotlist    = browser_tracks$data,
    rel_heights = track_heights,
    ncol        = 1,
    align       = "v",
    axis        = "rl"
  )
  
  final_plots
}



# ---- Figure 1 ----

# Function to create 5' metaplots
create_5_metaplots <- function(input_df, plot_title, stranded = F, 
                               plot_SEM = T, trans = 100, 
                               plot_levels = c("Control", "+TFIISDN"),
                               x_breaks = c(-2, 0, 2), 
                               x_labs = c("-2 kb", "TSS", "+2 kb"), ...) {
  
  # Simplify plot_title for matching
  match_str <- str_split(plot_title, " |,", simplify = T)[1, 1] %>%
    str_remove("-seq$") %>%
    str_remove("^m") %>%
    str_c("^", .)
  
  # For complicated titles input match string followed by title 
  plot_title <- plot_title %>%
    str_remove(str_c(match_str, ","))
  
  # Rename keys
  res <- input_df %>%
    filter(grepl(match_str, key)) %>%
    mutate(type = plot_title) %>%
    rename_key()
  
  # Identify plot levels
  if (all(grepl("sense", res$key))) {
    plot_levels <- c(
      "Control anti-sense",  "Control sense", 
      "+TFIISDN anti-sense", "+TFIISDN sense"
    )
  }
  
  # Create 5' metaplot
  res <- res %>%
    create_metaplots(
      stranded    = stranded,
      plot_SEM    = plot_SEM,
      plot_levels = plot_levels,
      trans       = trans,
      x_breaks    = x_breaks,
      x_labs      = x_labs,
      v_line      = 0,
      size        = 2,
      ...
    ) +
    facet_wrap(~type) +
    labs(y = "Mean relative signal") +
    theme(axis.text.x = element_text(hjust = c(0.5, 0.5, 0.7)))
  
  res
}



# ---- Figure 2 ----

# Function to create 3' metaplots
create_3_metaplots <- function(input_df, plot_title, plot_SEM = T, line_size = 3,
                               trans = 1000, plot_levels = c("Control", "+TFIISDN"), 
                               y_breaks = seq(0, 4, 0.6), ...) {
  
  # Simplify plot_title for matching
  match_str <- str_split(plot_title, " |,", simplify = T)[1, 1] %>%
    str_remove("-seq$") %>%
    str_remove("-P|^m") %>%
    str_c("^", .)
  
  # For complicated titles input match string followed by title 
  plot_title <- plot_title %>%
    str_remove(str_c(match_str, ","))
  
  # Create 3' metaplot
  res <- input_df %>%
    filter(grepl(match_str, key)) %>%
    mutate(type = plot_title) %>%
    rename_key() %>%
    create_metaplots(
      stranded    = F,
      plot_SEM    = plot_SEM,
      plot_levels = plot_levels,
      trans       = trans,
      x_breaks    = c(0, 5),
      x_labs      = c("pAS", "+5 kb"),
      y_breaks    = y_breaks,
      v_line      = 0,
      # n_label     = c(0.82, 0.8),
      n_label     = c(0.85, 0.7),
      size        = line_size,
      ...
    ) +
    facet_wrap(~type) +
    labs(y = "Mean signal (RPKM)") +
    theme(
      axis.title.y = element_text(size = 14, color = "black", face = "bold"),
      axis.text    = element_text(size = 10),
      axis.text.x  = element_text(hjust = c(0.5, 0.7))
      # axis.text.y  = element_text(size = 10, color = "black", face = "bold")
    )
  
  res
}



# ---- Figure 3 ----

# Function to identify pol II waves based on background signal
find_waves <- function(input_df, gene_list = NULL, max_dist = 60, sd_lim = 8) {
  
  # Print memory usage and system time
  system.time({
    mem_change({
      
      # Calculate gene lengths (Kb)
      if (!is.null(gene_list) & !"gene_len" %in% names(gene_list)) {
        gene_list <- gene_list %>%
          mutate(gene_len = round((end - start) / 1000), digits = 1) %>%
          dplyr::select(name, gene_len)
      }
      
      # Smooth data
      res <- input_df %>%
        ungroup()
      
      # Create nested data.frame
      if (!"data" %in% names(res)) {
        res <- res %>%
          nest(-key)
      }
      
      # Identify wave coordinates
      res <- res %>%
        mutate(data = map(data, ~ {
          
          max_wave_dist = max_dist - 5
          waves <- .x
          
          # Remove genes that are too short
          if (!is.null(gene_list)) {
            waves <- waves %>%
              left_join(gene_list, by = "name") %>%
              na.omit() %>%
              filter(gene_len >= max_dist)
          }
          
          # Remove windows that are too far downstream
          waves <- waves %>%
            filter(win_dist <= max_dist) %>%
            mutate(win_type = ifelse(win_dist <= max_wave_dist, 1, 2)) %>%                         
            group_by(name, win_type) %>%
            
            # Calculate mean count and sd for each timepoint
            mutate(
              mean_count = mean(count), 
              sd_count   = sd(count),
              limit      = mean_count + (sd_lim * sd_count)
            ) %>%
            group_by(name) %>% 
            mutate(limit = ifelse(win_dist <= max_wave_dist, min(limit), limit)) %>%  
            
            # Identify most distal bin where signal is above the limit
            filter(count > limit) %>%
            arrange(desc(win_dist)) %>%
            dplyr::slice(1) %>% 
            ungroup() %>%
            rename(wave_coord = win_dist)
          
          waves
        }))
      
      res <- res %>%
        unnest() %>%
        dplyr::select(  # Want to keep extra columns from input    
          -win_id,   -count, 
          -win_len,  -gene_len, 
          -win_type, -mean_count,     
          -sd_count, -limit
        )
      
    }) %>%
      cat(sep = "", "Memory usage:\n  ", ., " B\n\nSystem time:\n")
    
  }) %>%
    print()
  
  res
}

# Function to calculate transcription rates
calculate_rates <- function(input_df, rate_col = "m20_m10_rate") {
  
  # Function to calculate elongation rates
  calc_rates <- function(input_df, wave_col = "wave_coord") {
    
    # Function to extract timepoint from key
    extract_tm <- function(input_name) {
      res <- input_name %>%
        str_split("_") %>%
        unlist() %>%
        grep("m[0-9]", ., value = T)
      
      res
    }
    
    # Extract timepoint from key
    res <- input_df %>%
      mutate(DRB = extract_tm(key)) %>%
      dplyr::select(name, DRB, wave_col)
    
    # Create named list of timepoints
    tms <- unique(res$DRB)
    
    tm_list <- map(tms, ~ {
      tm <- str_extract_all(.x, "[0-9]") %>% 
        unlist() %>% 
        str_c(collapse = "") %>% 
        as.numeric()
      
      tm
    }) %>%
      set_names(tms)
    
    res <- res %>% 
      spread(DRB, wave_col)
    
    # Calculate rates for each combination of timepoints
    for (i in seq_along(tm_list)) {
      tm_1  <- tm_list[[i]]
      col_1 <- names(tm_list)[i] 
      
      other_tms <- tm_list[ tm_1 > tm_list ] %>%
        unlist()
      
      if (length(other_tms) > 0) {
        for (j in seq_along(other_tms)) {
          tm_2  <- other_tms[[j]]
          col_2 <- names(other_tms)[j]
          
          tm <- tm_1 - tm_2
          col_name <- str_c(col_1, col_2, "rate", sep = "_") %>%
            sym()
          
          res <- res %>%
            mutate(!!col_name := (!!sym(col_1) - !!sym(col_2)) / tm)
        }
      }
    }
    
    res  
  }
  
  # Function to filter each trxn rate column
  filt_rates <- function(input_df, rate_col) {
    
    rate_syms <- syms(rate_col)
    res       <- input_df
    
    # Filter each trxn rate column
    for (i in seq_along(rate_syms)) {
      filt_expr <- expr(!!rate_syms[[i]] > 0) 
      
      res <- res %>%
        filter(eval(filt_expr))
    }
    
    res
  }
  
  # Extract column names for timepoints
  tm_cols <- rate_col %>%
    str_split("_") %>%
    unlist() %>%
    grep("[0-9]", ., value = T)
  
  # Divide key into type and treatment
  res <- input_df %>%
    filter(grepl("m[0-9]", key)) %>%
    separate(
      col    = key, 
      sep    = "_", 
      into   = c("type", "treatment"), 
      extra  = "drop", 
      remove = F
    ) %>%
    group_by(type, treatment) %>%
    
    # Calculate rates
    do({calc_rates(.)}) %>%
    ungroup() %>%
    
    # Filter rates
    filt_rates(rate_col = rate_col) %>%
    dplyr::select(type, treatment, name, tm_cols, rate_col) %>%
    
    # Remove genes that are not present in both treatments
    group_by(type) %>%
    mutate(n_treat = n_distinct(treatment)) %>%
    group_by(type, name) %>%
    filter(n() == n_treat) %>%
    
    # Get gene count
    group_by(type, treatment) %>%
    mutate(n_gene = n_distinct(name)) %>%
    ungroup()
  
  # Print elongation rate stats
  tm_1 <- sym(tm_cols[1])
  tm_2 <- sym(tm_cols[2])
  rate_sym <- sym(rate_col)
  
  res %>%
    group_by(type) %>%
    summarize(
      n_gene = n(),
      !!tm_2     := median(!!tm_2),
      !!tm_1     := median(!!tm_1),
      !!rate_sym := median(!!rate_sym)
    ) %>%
    print()
  
  res
}

# Function to run DESeq2
run_DESeq2 <- function(input_df, alpha = 0.05, ...) {

  # Function to assign condition labels to DESeq design table 
  get_cond <- function(input_df, untreated, treated) {
    in_names <- names(input_df)
    
    res <- map(in_names, ~ {
      if (str_detect(.x, untreated)) {
        return("untreated")
      } else if (str_detect(.x, treated)) {
        return("treated")
      }
    }) %>%
      unlist()
    
    res
  }
  
  # Function to create DESeq design table
  make_design <- function(input_df) {
    sample_num <- input_df %>%
      colnames() %>% 
      length()
    
    res <- data.frame(
      row.names = names(input_df),
      condition = get_cond(input_df, ...),
      type      = c(rep("single-end", times = sample_num))
    )
    
    res
  }
  
  # Added rownames to input_df
  rownames(input_df) <- input_df$name
  
  input_df <- input_df %>% 
    dplyr::select(-name)
  
  # Create DESeq design table
  DESeq_design <- make_design(input_df)
  
  DESeq_object <- DESeqDataSetFromMatrix(
    countData = input_df,
    colData   = DESeq_design,
    design    = ~condition
  )
  
  # Run DESeq2
  DESeq_res <- DESeq(DESeq_object)
  
  # Set order of conditions
  DESeq_res$condition <- fct_relevel(
    DESeq_res$condition, 
    c("untreated", "treated")
  )
  
  # Get DESeq2 results
  DESeq_res <- DESeq_res %>%
    results(alpha = alpha)
  
  # Convert results to data.frame
  col_order <- c(
    names(genes), 
    names(DESeq_res)
  ) %>%
    grep("score", ., value = T, invert = T) %>%
    grep("symbol", ., value = T, invert = T)
  
  DESeq_df <- DESeq_res %>%
    data.frame() %>%
    mutate(name = rownames(DESeq_res)) %>%
    left_join(genes, by = "name") %>%
    mutate(gene_len = end - start) %>%
    dplyr::select(col_order, gene_len)
  
  DESeq_df
}

# Function to divide data.frame into groups 
get_groups <- function(input_df, num_groups = 4) {
  
  row_nums   <- seq_along(input_df$name)
  row_tot    <- length(row_nums)
  group_size <- floor(row_tot / num_groups)
  
  res <- input_df %>%
    mutate(key = NA)
  
  for (i in seq(1 : num_groups)) {
    min_row <- group_size * (i - 1) + 1
    max_row <- group_size * i
    
    group_coords <- row_nums[min_row : max_row]
    
    if (i == num_groups & row_tot > max_row) {
      extra_rows   <- row_nums[max_row + 1 : row_tot]
      group_coords <- c(group_coords, extra_rows)
    }
    
    res <- res %>%
      mutate(key = ifelse(row_nums %in% group_coords, as.character(i), key))
  }
  
  res
}

# Funciton to calculate p-values for metaplots
calc_pvals <- function(input_df, rep_mean = F, 
                       key_names = c("key", "rep", "strand")) {
  
  keys <- input_df$key %>%
    unique()
  
  if (rep_mean) {
    no_rep_pos   <- grep("rep", key_names, invert = T)
    no_rep_names <- key_names[no_rep_pos]
    split_keys   <- str_split(keys, "_")
    
    no_rep_keys <- map(split_keys, ~ {
      .x[no_rep_pos] %>%
        str_c(collapse = "_")
    }) %>%
      unlist() %>%
      unique()
    
    data_1 <- sym(no_rep_keys[1])
    data_2 <- sym(no_rep_keys[2])
    
    res <- input_df %>%
      separate(key, sep = "_", into = key_names) %>%
      unite("key", no_rep_names, sep = "_") %>%
      spread(key, count) %>%
      group_by(rep, win_id) %>%
      
      # Calculate p-values
      summarize(
        p_val = tidy(t.test(!!data_1, !!data_2))$p.value
      ) %>%
      ungroup() %>%
      mutate(
        p_val = -log10(p_val),
        n_rep = n_distinct(rep)
      ) %>%
      
      # Calculate mean and SEM for replicates
      group_by(win_id, n_rep) %>%
      summarize(
        SEM   = sd(p_val),
        p_val = mean(p_val)
      ) %>%
      ungroup() %>%
      mutate(SEM = SEM / sqrt(n_rep))
    
  } else {
    data_1 <- sym(keys[1])
    data_2 <- sym(keys[2])
    
    res <- input_df %>%
      spread(key, count) %>%
      group_by(win_id) %>%
      summarize(
        p_val = tidy(t.test(!!data_1, !!data_2))$p.value
      ) %>%
      mutate(p_val = -log10(p_val))
  }
  
  res
}



# ---- Figures 4 and 5 ----

# Function to merge and normalize pause density tables 
norm_pause_density <- function(input_df, length_norm = T, ...) {

  res <- input_df
  
  # Normalize pause counts by window length
  if (length_norm) {
    res <- res %>%
      mutate(data = map2(data, key, ~{
        if (grepl("pauses", .y))
          mutate(.x, count = count / (end - start))
        else .x
      }))
  }
    
  # Merge data.frames
  res <- res %>%
    merge_wins(...) %>%
    unnest() %>%
  
    # Calculate mean NET-seq signal for each window
    separate(key, sep = "_", into = c("type", "key", "rep", "region")) %>%
    group_by(type, key, name, win_id, region) %>%
    summarize(count = mean(count)) %>%
    ungroup() %>%
    spread(type, count) %>%
    na.omit() %>%

    # Add pseudo count
    mutate(zero = ifelse(NET == 0, T, F)) %>%
    group_by(key, name, zero) %>%
    mutate(min_val = min(NET)) %>%
    ungroup() %>%
    
    # Normalize pause density by NET-seq signal
    group_by(key, name) %>%
    mutate(
      NET  = ifelse(NET == 0, max(min_val) / 2, NET),
      norm = pauses / NET
    ) %>%
    ungroup() %>%
    dplyr::select(-zero, -min_val) %>%
    gather(type, count, NET, pauses) %>%
    unite(key, type, key, region, sep = "_")
  
  res
}

# Function to plot mean pause density
create_density_metaplots <- function(input_df, y_title, show_legend = F, ...) {
  
  # Function to plot pause density
  plot_density <- function(input_df, gene_list, ...) {
    res <- input_df %>%
      semi_join(gene_list, by = "name") %>%
      calc_mean_signal(
        data_col   = "norm",
        group_cols = c("type", "key", "region")
      ) %>%
      rename_key() %>%
      create_metaplots(
        data_col    = "norm",
        v_line      = 0,
        size        = 3,
        ...
      ) +
      theme(axis.text.x     = element_text(hjust = c(0.5, 0.7))) +
      labs(y = y_title)
    
    res
  }
  
  # Create 5' pause density plots
  # Nonoverlapping genes >5 kb long
  genes_5 <- gene_lists[["sep_0k"]] %>%
    filter(end - start > 5000) %>%
    dplyr::select(name)
  
  den_5 <- input_df %>%
    filter(grepl("pauses_.+_5", key)) %>%
    plot_density(
      gene_list   = genes_5,
      x_breaks    = c(0, 5),
      x_labs      = c("TSS", "+5 kb"),
      ...
    ) +
    theme(legend.position = "none")
  
  # Create 3' pause density plots
  # Genes >5 kb separated and >2 kb long
  genes_3 <- gene_lists[["sep_5k"]] %>%
    filter(end - start > 2000) %>%
    dplyr::select(name)
  
  den_3 <- input_df %>%
    filter(grepl("pauses_.+_3", key)) %>%
    plot_density(
      gene_list   = genes_3,
      x_breaks    = c(-2, 0, 5),
      x_labs      = c("-2 kb", "pAS", "+5 kb"),
      ...
    ) +
    theme(
      axis.title   = element_blank(),
      axis.text.x  = element_text(hjust = c(0.5, 0.5, 0.7)),
      axis.line.y  = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  if (!show_legend) {
    den_3 <- den_3 +
      theme(legend.position = "none")
  }
  
  res <- plot_grid(
    den_5, den_3,
    nrow  = 1,
    align = "h"
  )
  
  res
}

# Function to create metaplots for region around pause sites
create_pause_metaplots <- function(input_df, ...) {
  res <- input_df %>%
    rename_key() %>%
    create_metaplots(
      plot_SEM    = T,
      plot_levels = c("Control", "+TFIISDN"),
      trans       = 1000,
      x_lim       = c(-0.016, 0.015),
      x_breaks    = seq(-0.01, 0.01, 0.01),
      x_labs      = c("-10 bp", "Pause", "+10 bp"),
      v_line      = 0,
      n_label     = c(0.82, 0.8),
      ...
    ) + 
    labs(y = "Mean signal (RPKM)")
  
  res
}

# Function to create sequence logos
create_logo <- function(input_freq, plot_colors, x_lim = NULL, y_lim = NULL,
                        x_breaks = NULL, x_labs = NULL, v_line = NULL) {
  
  # Set color scheme
  plot_cols <- make_col_scheme(
    chars = c("A", "T", "G", "C"), 
    cols  = plot_colors
  )
  
  # Simplify table 
  freq_df <- input_freq %>% 
    mutate(
      A    = round(A, 3), 
      `T`  = round(`T`, 3), 
      G    = round(G, 3), 
      C    = round(C, 3),
      tot  = (A + `T` + G + C),
      diff = 1 - tot,
      C    = ifelse(tot != 1, C + diff, C)
    ) %>%
    dplyr::select(-tot, -diff) %>%
    gather(nuc, value, -position)
  
  # Set x-limits 
  if (!is.null(x_lim)) {
    x_min <- x_lim[1]
    x_max <- x_lim[2]
    
    freq_df <- freq_df %>% 
      filter(position >= x_min & position <= x_max)
  }
  
  freq_df <- freq_df %>% 
    spread(position, value) 
  
  # Add rownames
  rownames(freq_df) <- freq_df$nuc
  
  freq_df <- freq_df %>% 
    dplyr::select(-nuc)
  
  # Create sequence logo
  res <- ggseqlogo(as.matrix(freq_df), col_scheme = plot_cols) +
    scale_y_continuous(breaks = seq(0, 2, 0.5))
  
  # Change x-axis labels
  if (!is.null(x_breaks)) {
    x_breaks <- seq(1 : (x_lim[2] - x_lim[1] + 1))
    
    res <- res +
      scale_x_continuous(
        breaks = x_breaks, 
        labels = x_labs
      )
  }
  
  # Add vertical dotted line
  if (!is.null(v_line)) {
    res <- res + 
      geom_vline(
        xintercept = v_line, 
        linetype   = 2, 
        size       = 1
      )
  }
  
  # Set y-limits
  if (!is.null(y_lim)) {
    res <- res + 
      coord_cartesian(ylim = y_lim)
  }
  
  res + theme_info
}

# Function to export pause bed file
write_pause_bed <- function(input_df, file_name) {
  input_df %>%
    mutate(name = str_remove(name, "\\*[:digit:]+,[:digit:]+$")) %>%
    left_join(pause_lists[["H2O_sep"]], by = "name") %>%
    dplyr::select(chrom, start, end, name, score = 7, strand) %>%
    bed_sort() %>%
    unique() %>%
    write_tsv(path = str_c(R_list_dir, file_name))
}

# Function to run cross correlation analysis 
calc_cross_corr <- function(input_df, gene_list) {
  
  # Function to create vector of read counts for every base pair position on every gene 
  getCounts <- cppFunction(
    'std::vector<double> getCounts(DataFrame inputBed) {
    
      std::vector<int> startCoords = inputBed["start"];
      std::vector<int> endCoords   = inputBed["end"];
      std::vector<double> counts   = inputBed["count"];
    
      int N = startCoords.size();
    
      std::vector<double> countVec;
    
      for (int i = 0; i < N; i++) {
        int coordSize = endCoords[i] - startCoords[i];
    
        for (int j = 0; j < coordSize; j++) {
          countVec.push_back(counts[i]);
        }
      }
    
      return countVec;
    }'
  )
  
  # Create nested data.frame containing counts for each bp position
  res <- input_df %>%
    mutate(data = map(data, ~{
      .x %>%
        semi_join(gene_list, by = "name") %>%
        group_by(name) %>%
        filter(sum(count) > 0) %>%
        nest() %>%
        mutate(data = map(data, getCounts)) %>%
        ungroup()
    }))
  
  # Remove genes that are not present in both datasets 
  res <- res$data %>%
    purrr::reduce(left_join, by = "name") %>%
    mutate(data_type = map2(data.x, data.y, ~{
      str_c(typeof(.x), typeof(.y))
    })) %>%
    unnest(data_type) %>%
    filter(!grepl("NULL", data_type)) %>%
    dplyr::select(-data_type) %>% 
    
    # Perform cross correlation analysis 
    mutate(corr = map2(data.x, data.y, ~{
      ccf(.x, .y, plot = F, type = "correlation") %>%
        tidy() %>%
        filter(lag >= -20 & lag <= 20)
    })) %>%
    dplyr::select(name, corr) %>%
    unnest()
  
  res
}



# ---- Figure 6 ----

# Function to identify refernce Cq values
identify_ref_Cq <- function(input_df) {
  
  # Calculate mean Cq value across tech reps for each individual qPCR run
  res <- input_df %>%
    ungroup() %>%
    separate(target, sep = "_", into = "target", extra = "drop") %>%
    group_by(target, date, rep, sample) %>%
    summarize(Cq = mean(Cq)) %>%
    ungroup() %>%
    
    # Identify refence Cq values
    spread(target, Cq) %>% 
    gather(target, Cq, -sample, -date, -rep, -ACTB) %>%
    dplyr::rename(ref = ACTB) %>%
    na.omit() %>% 
    unite(Cq, Cq, ref, sep = ",")
  
  res
}

# Function to calculate delta delta Cq
calc_ddCq <- function(input_df) {
  
  # Create columns for sample and control, target and reference Cq values
  res <- input_df %>%
    dplyr::select(target, date, rep, sample, Cq, con) %>%
    separate(Cq, sep = ",", into = c("sam_targ", "sam_ref"), convert = T) %>%
    separate(con, sep = ",", into = c("con_targ", "con_ref"), convert = T) %>%
    
    # Calculate delta-delta-Cq values    
    mutate(
      "sam_targ-ref"        = sam_targ - sam_ref,
      "con_targ-ref"        = con_targ - con_ref, 
      "delta_sam-delta_con" = `sam_targ-ref` - `con_targ-ref`,
      fold_change           = (2 ^ (-`delta_sam-delta_con`))
    ) %>%
    na.omit() %>% 
    
    # Calculate mean fold change across qPCR runs (tech reps) for each biological replicate
    group_by(target, sample, rep) %>%
    summarize(fold_change = mean(fold_change)) %>%
    
    # Calculate mean fold change and SEM for biological replicates
    group_by(target, sample) %>% 
    summarize(
      n_rep   = n_distinct(rep),
      mean_FC = mean(fold_change),
      SEM     = sd(fold_change) / sqrt(n_rep)
    ) %>%
    ungroup() %>%
    
    # Remove targets that only have one biological replicate
    na.omit() %>%
    separate(sample, sep = "_", into = c("key", "tm"))
  
  res
}

# Function to calculate qPCR fold change 
calc_qPCR_FC <- function(input_df) {
  
  # Identify reference Cq values
  res <- input_df %>%
    identify_ref_Cq() %>%
    
    # Identify -DMOG Cq values
    separate(sample, sep = "_", into = c("key", "tm")) %>% 
    spread(tm, Cq) %>%
    mutate(con = Unt) %>% 
    gather(tm, Cq, -target, -key, -rep, -date, -con) %>%
    unite(sample, key, tm, sep = "_") %>%
    
    # Calculate delta delta Cq
    calc_ddCq()
  
  res
}  

# Function to create bargraphs for qPCR data
create_qPCR_bargraphs <- function(input_df, gene_targets, plot_colors, tm_levels) {
  
  # Generate plot labels
  tm_labs <- tm_levels %>%
    str_replace("30m", "0.5h") %>%
    str_replace("h", " h")
  
  key_labs <- c("Control", "+TFIISDN") %>%
    map(add_subscript)
  
  # Set plot levels
  res <- input_df %>%
    filter(
      target %in% gene_targets,
      tm     %in% tm_levels
    ) %>%
    mutate(
      target = fct_relevel(target, gene_targets),
      key    = fct_relevel(key, c("H2O", "Dox")),
      tm     = fct_relevel(tm, tm_levels)
    ) %>%
    
    # Create bargraphs with error bars
    ggplot(aes(tm, mean_FC, 
      ymin  = mean_FC, 
      ymax  = mean_FC + SEM, 
      color = key, 
      fill  = key
    )) +
    geom_bar(
      stat        = "identity", 
      width       = 0.8, 
      position    = position_dodge(width = 0.9), 
      show.legend = F
    ) +
    geom_errorbar(
      width    = 0.15, 
      position = position_dodge(width = 0.9)
    ) +
    
    # Set plot plot colors and labels
    scale_x_discrete(breaks = tm_levels, labels = tm_labs) +
    scale_colour_manual(
      values = plot_colors,
      labels = key_labs,
      guide  = guide_legend(override.aes = list(size = 3))
    ) +
    scale_fill_manual(
      values = plot_colors,
      labels = key_labs
    ) +
    
    # Adjust theme attributes
    theme_info +
    theme(
      strip.background = element_blank(),
      strip.text.x     = element_text(size = 18, face = "bold"),
      axis.title.x     = element_blank(),
      legend.text      = element_text(size = 16, face = "bold"),
      legend.position  = "right"
    ) +
    facet_wrap(~target, nrow = 2, scales = "free_y")
   
  res
}


