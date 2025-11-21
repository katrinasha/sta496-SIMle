page <- read_html('hkepd_dat.htm')
page %>% html_table(fill =TRUE)
html_nodes(page, 'table')
cat(as.character(page)[1:5000])
page %>% html_node('pre') %>% html_text()

library(readr)

# Path to your local file
f <- "hkepd_dat.htm"   # change path if needed

# 1) Read all lines with the correct encoding
lines <- readLines(f, encoding = "windows-1252")

# 2) Keep only the lines inside the <PRE> ... </PRE> block
start <- grep("<pre>",  lines, ignore.case = TRUE)[1]
end   <- grep("</pre>", lines, ignore.case = TRUE)[1]
pre_lines <- lines[(start + 1):(end - 1)]

# 3) Find the first data row (i.e., the first line that starts with a number)
first_data <- which(grepl("^\\s*[-+]?\\d", pre_lines))[1]

# 4) Collapse the data rows to text and read as whitespace-delimited table
dat_text <- paste(pre_lines[first_data:length(pre_lines)], collapse = "\n")

hk <- read.table(text = dat_text, header = FALSE, sep = "", strip.white = TRUE,
                 stringsAsFactors = FALSE)

# 5) Add column names
names(hk) <- c("SO2","NO2","RSP","O3","Temp","Humidity",
               "Circulatory","Respiratory","Total","w1","w2","w3","w4","w5","w6")

# Quick sanity checks
str(hk)
head(hk)
write.csv(hk, 'hkepd_dat.csv')