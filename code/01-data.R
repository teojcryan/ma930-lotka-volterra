### load data ----------------
lynx_hare_df = data.frame(
  read.csv('data/hare_lynx_data.csv',
           fileEncoding="UTF-8-BOM")) %>%
  rename(prey = 2, pred = 3)

data = lynx_hare_df[1:68,]
times = data$year