library(dplyr)
library(tidyr)
library(ggplot2)
library(odbc)


query <- "SELECT
    Accession,
    protein_name,
    ec_numbers,
    inference_level,
    organism,
    infoplace
    FROM
    annotation_table"


con <- DBI::dbConnect(odbc(), # connects us to the management system that we want to analyse our data
                      Driver = "SQL Server", #Type of driver
                      Server = "DDI01",
                      Database = "HCP_Annotator_DB", #Database name
                      Trusted_Connection = "True")



result <- dbGetQuery(con, query)

#dbDisconnect(con)

add_protein_class <- function(df, critical_class_label = "I") {
  df$Class <- NA
  df$Reasoning <- NA
  
  # Funktion til at bestemme om EC-nummeret matcher kriterierne for the critical class
  
  is_critical_enzyme <- function(ec_number) {
    if (!is.na(ec_number)) {
      return(startsWith(ec_number, "3.2.") | startsWith(ec_number, "3.4.") | startsWith(ec_number, "3.1.1."))
    } else {
      return(FALSE)
    }
  }    
  
  # 1. Separate rows based on comma-separated EC numbers
  df_separated <- df %>%
    separate_rows(ec_numbers, sep = ",\\s*")
  
  
 
  # 2. Determine the individual class for each EC number
  df_classified <- df_separated %>%
    rowwise() %>%
    mutate(
      Individual_Class = ifelse(is_critical_enzyme(ec_numbers), critical_class_label, "III")
    )

  
  # 3. Determine the final class and reasoning per Accession
  df_final <- df_classified %>%
    group_by(Accession) %>%
    summarise(
      Has_Critical = any(Individual_Class %in% c(critical_class_label)),
      Class = ifelse(Has_Critical, critical_class_label, "III"),
      Reasoning = ifelse(Has_Critical, "Critical Enzyme Activity", "No critical enzymatic activity detected"),
      ec_numbers = paste(ec_numbers, collapse = ", "),
      # Keep other relevant columns (you might need to adjust this based on your actual DF)
      across(everything(), ~ first(.)),
      .groups = "keep"
    ) %>%
    select(-Has_Critical,-Individual_Class)
  
  return(df_final)
}





query <- "SELECT
    Accession,
    Positive_High_assay_number,
    Positive_Intermediate_assay_number,
    Positive_assay_number,
    Positive_Low_assay_number,
    Negative_assay_number
    FROM
    iedb_table"

immune_table <-  dbGetQuery(con, query) %>%
  mutate(Allpositive =  Positive_High_assay_number + Positive_Intermediate_assay_number + Positive_assay_number + Positive_Low_assay_number ) %>%
  filter(Allpositive > Negative_assay_number)




immune_prot <- result[result$Accession %in% immune_table$Accession, ,drop = F]
immune_prot$Class <- "II"
immune_prot$Reasoning <- "Immunogenic Activity in ImmunoAssays"






protein_level <- all_prot_info_update_1805[all_prot_info_update_1805$inference_level == "1: Evidence at protein level", ,drop = F]
protein_level <- add_protein_class(protein_level, critical_class_label = "I")

other_levels <- all_prot_info_update_1805[all_prot_info_update_1805$inference_level != "1: Evidence at protein level", ,drop = F]
other_levels_t <- add_protein_class(other_levels, critical_class_label = "II")

all_prot <- rbind(protein_level, other_levels_t)


all_prot <- all_prot_info_update_1805
immune_prot <- HCPAnnotatorDB[["iedb_table"]]

all_prot[all_prot$Accession %in% immune_prot$Accession, "Class"] <- "II"

all_prot[all_prot$Accession %in% immune_prot$Accession & all_prot$Reasoning %in% "No critical enzymatic activity detected", "Reasoning"] <- "Immunogenic Potential"

all_prot[all_prot$Accession %in% immune_prot$Accession & all_prot$Reasoning %in% "Critical Enzyme Activity", "Reasoning"] <- "Critical Enzyme Activity, Immunogenic Potential"

all_prot[all_prot$Reasoning %in% "No critical enzymatic activity detected", "Reasoning"] <- "No direct risk Observed"

all_prot_info_update_1805 <- all_prot

library(dplyr)


table(all_prot$Class)
table(all_prot$Reasoning)

sum(grepl("Immunogenic Activity in ImmunoAssays",all_prot_info_update_1805_info_update_1805$Reasoning))
table(all_prot_info_update_1805_info_update_1805$Reasoning)

all_prot[all_prot$Accession %in% c("A0A9J7GF01", "A0A9J7GF01","A0A8C2QK35","G3HR95","G3I4W7","G3INC5","A0A3L7IKX6","A0A9J7J6W2","G3H6V7","A0A061IKA1"),"Class"] <- "I"
all_prot[all_prot$Accession %in% c("A0A9J7GF01", "A0A9J7GF01","A0A8C2QK35","G3HR95","G3I4W7","G3INC5","A0A3L7IKX6","A0A9J7J6W2","G3H6V7","A0A061IKA1"),"Reasoning"] <- "Critical Enzymatic Activity (Literature)"


all_prot[all_prot$Accession %in% c("A0A061I4M6","G3I6T1"),"Class"] <- "I"
all_prot[all_prot$Accession %in% c("A0A061I4M6","G3I6T1"),"Reasoning"] <-"Immunogenic Potential (Clinical Trial)"


  

all_prot_info <- merge(all_prot, df_potential_impact_manual,  by.x = "Accession",  by.y = "V1", all.x = T)


for_just_plot <- merge(all_prot, result %>% select(Accession,organism,infoplace), by = "Accession", all.x=T)

#setwd("~/Master/Tasks/Scripts_RObjects_January")

#save(all_prot, file = "All_new_columns_updated.RData")



#for_just_plot <- merge(test_1, result %>% select(Accession, organism), by = "Accession", all.x = T)
for_just_plot_data <- for_just_plot[for_just_plot$infoplace != "" &  for_just_plot$infoplace != "ELISA-MS" & for_just_plot$infoplace != "Not Available" ,, drop = F]
for_just_plot_data <- for_just_plot_data [for_just_plot_data$infoplace != "Not Available",,drop = F]

plot_data_1 <- for_just_plot_data %>%
  separate_rows(Reasoning, sep = ",\\s*") %>%
  group_by(organism,Class, Reasoning) %>%
  summarise(total_groups = n())




plot_data_1 <- plot_data_1 %>%
  group_by(organism) %>%
  mutate(total = sum(total_groups),percentage = round(total_groups / total * 100, 2))

table(all_prot$Class)


plot_data_1$organism <- factor(plot_data_1$organism, levels = c("CHO",
                                                                "Human",
                                                                "Mesocricetus auratus",
                                                                "Chlorocebus sabaeus",
                                                                "Gallus gallus",
                                                                "Coturnix japonica",
                                                                "Nicotiana tabacum" ,
                                                                "Physcomitrium patens",
                                                                "Drosophila melanogaster",
                                                                "Spodoptera frugiperda",
                                                                "Acinetobacter baumannii",
                                                                "Shigella flexneri",
                                                                "Staphylococcus aureus",
                                                                "Pseudomonas aeruginosa",
                                                                "Streptoverticillium mobaraense",
                                                                "Escherichia coli",
                                                                "Saccharomyces cerevisiae"
))

plot_data_1 <- plot_data_1[plot_data_1$organism != "Shigella flexneri",, drop =F]


border_color <- "#4A5568"  # Darker gray for outlines
colors_p <- c("#0073b7", "#cdd4dd", "#f39c12","#9baabb", "#FFC945","#ffffff", "#375577","#b27300","#6C7E8C", "#062b56")


table(all_prot$Class)
table(all_prot$Reasoning)

plot_data_1 <- plot_data_1[!is.na(plot_data_1$Class),, drop = F]
plot1 <- ggplot(plot_data_1, aes(x = organism, y = percentage, fill = Class)) +
  geom_bar(stat = "identity", color = border_color, position = "stack") +
  scale_fill_manual(values = c("#FFC945","#cdd4dd","#455770"))+
  
  geom_text(
    aes(label = ifelse(percentage >= 4 , 
                       percentage,
                       "")), 
    position = position_stack(vjust = 0.5),
    size = 5,
    color = "black"
  ) +
  labs(title = "Organism Enzymatic Number", fill = "EC Numbers") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 16))

# 2. Create the second plot: Stacked bar plot of EC number distribution per organism (filtered)
# 2.1 Filter out rows where ec_numbers is empty
filtered_annotation_table <- annotation_table %>%
  filter(!is.na(ec_numbers) & ec_numbers != "")

# 2.2 Extract the first EC number
filtered_annotation_table <- filtered_annotation_table %>%
  mutate(first_ec = gsub("\\..*", "", ec_numbers))

# 2.3 Summarize the data: count occurrences of each `first_ec` per `organism`
summary_data_2 <- filtered_annotation_table %>%
  group_by(organism, first_ec) %>%
  filter(!is.na(first_ec) & first_ec != ""& first_ec != "Not Available") %>%
  dplyr::count() %>%
  ungroup() %>%
  # 2.4 Calculate percentage per organism
  group_by(organism) %>%
  mutate(percentage = n / sum(n)*100) %>%
  ungroup()



############################ Class vs Reasoning PLot
plot_data_1 <- for_just_plot_data %>%
  separate_rows(Reasoning, sep = ",\\s*") %>%
  group_by(organism, Class, Reasoning) %>%
  summarise(count = n(), .groups = "drop")

plot_data_for_viz <- plot_data_1 %>%
  group_by(Class) %>%
  mutate(total_in_class = sum(count)) %>%
  ungroup() %>%
  group_by(Class, Reasoning) %>%
  summarise(
    sum_count_per_reasoning = sum(count),
    total_in_class = first(total_in_class),
    .groups = "drop"
  ) %>%
  mutate(percentage = sum_count_per_reasoning / total_in_class) %>%
  arrange(Class, Reasoning) 



sample_data_for_plot <- plot_data_1 %>%
  # --- Debugging Step 3: Check 'count' just before the error-causing line ---
  # If you still suspect 'count' is the issue, add:
  # { if (!is.numeric(.$count)) stop("Error: 'count' column is not numeric before sum!") } %>%
  group_by(Class) %>%
  mutate(total_in_class = sum(count)) %>% # This is where the error occurs
  ungroup() %>%
  group_by(Class, Reasoning) %>%
  summarise(
    sum_count_per_reasoning = sum(count),
    total_in_class = first(total_in_class),
    .groups = "drop"
  ) %>%
  mutate(percentage = round(sum_count_per_reasoning / total_in_class * 100,1) ) %>%
  arrange(Class, Reasoning) # For consistent display order


sample_data_for_plot <- sample_data_for_plot[!(sample_data_for_plot$Class %in% "III"),]

# --- Plotting Function ---


border_color <- "#4A5568"  # Darker gray for outlines
colors_p <- c( "#cdd4dd","#0073b7", "#f39c12","#9baabb", "#ffffff", "#375577","#b27300","#6C7E8C", "#062b56")



plot_pie_charts_by_class <- function(data, title, threshold = 20) {
  p <- ggplot(data, aes(x = "", y = percentage, fill = Reasoning)) + # x = "" for pie chart trick
    geom_bar(stat = "identity", width = 1, color = border_color) + # width = 1 for full circle
    coord_polar("y", start = 0) + # Key for pie chart: polar coordinates
    labs(
      title = title,
      fill = "Reasoning Category" # Legend title
    ) +
    facet_wrap(~ Class, ncol = 3) + # Separate pie chart for each Class
    theme_void() + # Minimal theme for pie charts
    theme(
      axis.text = element_blank(), # No axis text for pie charts
      axis.title = element_blank(), # No axis titles
      panel.grid = element_blank(), # No grid lines
      legend.title = element_text(size = 14), # Legend title size
      legend.text = element_text(size = 12),  # Legend item text size
      plot.title = element_text(size = 18, hjust = 0.5), # Centered plot title
      strip.text = element_text(size = 14, face = "bold") # Facet label text size
    ) +
    # Conditional text labels for percentages on the pie slices
    geom_text(
      aes(label = ifelse(percentage >= threshold, percentage, "")),
      position = position_stack(vjust = 0.5), # Center text within slices
      size = 4, # Adjust text size
      color = "black"
    ) +
    scale_fill_manual(values=colors_p ) # Color palette for reasoning categories
  
  return(p)
}


sample_data_for_plot$Reasoning <- factor(sample_data_for_plot$Reasoning, levels =unique(sample_data_for_plot$Reasoning))

# --- Generate the plot using the sample data ---
pie_charts_plot <- plot_pie_charts_by_class(
  sample_data_for_plot,
  "Distribution of Reasoning within Protein Classes",
  threshold = 10 # Only label percentages 2% and above
)

print(grouped_plot)
