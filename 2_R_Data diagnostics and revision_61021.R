#################################################################
#################################################################
############# Wildfire smoke and Covid-19 case fatality #########
############# ratios in the San Francisco Bay Area: #############
############# a synthetic control analysis ######################
############# Version: June 10th, 2021 ########################## 
############# Code 2/3: Data Diagnostics and revision ########### 
################################################################# 
#################################################################



rm(list =ls())
gc()
library(foreign)
library(tidyverse)
library(data.table)
library(zoo)
library(psych)
library(lubridate)
library(gridExtra)
library(grid)
library(ggplot2)
options(scipen=999)
options(digits=5)

#setwd("C:/Users/annak/Dropbox/Projects/2020_San Diego/Wildfire smoke and COVID-19/")

setwd("F:/Lara/Projects/Wildfire smoke/")


##########################################################################################
### Original data
##########################################################################################

### Daily data
df1 <- read.dta("./Data/Diagnostics/Covid_smoke_by_county_diagnostics.dta") %>% 
  mutate(eligible = ifelse( (countyname=="Alameda County"         |
                               countyname=="Contra Costa County"  | 
                               countyname=="Marin County"         |
                               countyname=="San Francisco County" |
                               countyname=="Santa Clara County"   |
                               countyname=="Sonoma County") &   state=="CA", "Yes", "No")) %>% 
  mutate(week = lubridate::isoweek(date2))                       #add week number

#%>%         
#  mutate(cases = ifelse(countyname == "Marin County" & date == 22126, 160, cases)) %>% 
#  mutate(cases = ifelse(countyname == "Marin County" & date == 22098, 118, cases)) 
#
#df1_Marin <- df1 %>% 
#  filter(countyname == "Marin County")


fig <- df1 %>% filter(eligible == "Yes") %>% 
  filter(countyname == "Marin County") %>% 
  filter(week >= 20) %>% 
  ggplot(aes(x=date2, y=total_cases))+
  theme(axis.title.x=element_blank())+
  facet_wrap(~ countyname,  ncol = 6) + 
  geom_line()


### Total COVID deaths by county
df1_tot_covid <- df1 %>%                                       
  group_by(countyname, state) %>%
  summarise(tot_covid = sum(death)) %>% 
  ungroup()

### Weekly data
df1_weekly <- df1 %>%                                  #generate weekly data
  group_by(countyname, state, week) %>%             
  summarise(CFR_same_day = mean(CFR_same_day),         #calculate weekly mean CFR and smoke
            CFR_14days   = mean(CFR_14dys),
            CFR_7_14days = mean(CFR_7_14days),
            death        = mean(death),
            heavy_pct    = mean(heavy_pct),
            mobility1 = mean(populationstayingathome),
            mobility2 = mean(populationnotstayingathome)
            ) %>% 
  ungroup() %>% 
  left_join(df1_tot_covid) %>% 
  filter(tot_covid >= 100) %>%                          #filter to counties with over 100 total COVID deaths
  filter(week >=20) %>%                                 #restrict to week 20 and after
  mutate(time = week) %>%                           
  mutate(smoke_bin = ifelse(heavy_pct >= 70, 1, 0)) %>%  #generate variable for average weekly smoke coverage>70% 
  mutate(id = group_indices(., countyname, state))


df1_first_smoke <- df1_weekly %>%                        #find the first week with smoke>70% 
  group_by(id, smoke_bin) %>% 
  filter(time == min(time)) %>% 
  slice(1) %>%
  ungroup() %>% 
  filter(smoke_bin==1) %>% 
  mutate(first_smoke=1)

df1_weekly <- df1_weekly %>% 
  left_join(df1_first_smoke) %>% 
  mutate(smoke_week1 = ifelse(first_smoke == 1, time, NA)) %>% 
  mutate(eligible = ifelse( (countyname=="Alameda County"         |
                               countyname=="Contra Costa County"  | 
                               countyname=="Marin County"         |
                               countyname=="San Francisco County" |
                               countyname=="Santa Clara County"   |
                               countyname=="Sonoma County") &   state=="CA", "Yes", "No"))


### Daily data (Lara added)
df1_daily <- df1 %>%                                  
  left_join(df1_tot_covid) %>% 
  filter(tot_covid >= 100) %>%                          #filter to counties with over 100 total COVID deaths
  filter(week >=20) %>%                                 #restrict to week 20 and after
  mutate(time = date-21914) %>%                           
  mutate(smoke_bin = ifelse(heavy_pct >= 70, 1, 0)) %>% #generate variable for average weekly smoke coverage>70% 
  mutate(id = group_indices(., countyname, state))


df1_first_smoke <- df1_daily %>%                        #find the first week with smoke>70% 
  group_by(id, smoke_bin) %>% 
  filter(time == min(time)) %>% 
  slice(1) %>%
  ungroup() %>% 
  filter(smoke_bin==1) %>% 
  mutate(first_smoke=1)

df1_daily <- df1_daily %>% 
  left_join(df1_first_smoke) %>% 
  mutate(smoke_week1 = ifelse(first_smoke == 1, time, NA)) 


#########################################################################
### Generate smoothed data
#########################################################################

### 1. Generate daily cases and deaths
# 1.1  Daily cases and deaths are the change in total cases/deaths from the day before (same as Lara)
df2 <- df1 %>% 
  select_if(names(.)  %in%  c("countyname", "state", "day", "countyfips", "statefips", "date", "date2", "total_cases", "total_d", "light", "medium", "heavy", "light_pct", "medium_pct", "heavy_pct", "eligible", "week", "populationstayingathome", "populationnotstayingathome")) %>%   
  group_by(countyname, state) %>%
  arrange(day) %>%
  mutate(cases = total_cases - lag(total_cases, default = first(total_cases))) %>% 
  mutate(death = total_d - lag(total_d, default = first(total_d))) %>% 
  ungroup()

# 1.2 Identify negative daily cases and deaths
cases.neg <- df2 %>% 
  filter(cases<0) 
write.csv(cases.neg, "./Data/Diagnostics/data_negative_daily_cases.csv")

death.neg <- df2 %>% 
  filter(death<0) 
write.csv(death.neg, "./Data/Diagnostics/data_negative_daily_death.csv")

# 1.3 For now replace the negative daily entries with the values from previous day (instead of converting to missing). Note that some remain negative (we need to address this better)
df2 <- df2 %>% 
  group_by(countyname, state) %>%
  arrange(day) %>%
  mutate(cases = ifelse(cases<0, lag(cases, n=1), cases)) %>% #replace negative daily cases with values from day before
  mutate(death = ifelse(death<0, lag(death, n=1), death)) %>% #replace negative daily death with values from day before
  ungroup() 

# Drop counties with cases/deaths that are still negative
cases.neg <- df2 %>% 
  filter(cases<0) %>% 
  select(countyname, state) %>% 
  distinct %>% 
  mutate(miss.cases = 1)

death.neg <- df2 %>% 
  filter(death<0) %>% 
  select(countyname, state) %>% 
  distinct %>% 
  mutate(miss.death = 1)

df2 <- df2 %>% 
  left_join(cases.neg) %>% 
  left_join(death.neg) %>% 
  filter(is.na(miss.cases)) %>% 
  filter(is.na(miss.death)) %>% 
  select(!c(miss.cases, miss.death))

### 2. Smoooth the daily trends with a rolling sum and generate CFR measures
df2 <- df2 %>% 
  group_by(countyname, state) %>%      
  arrange(day) %>%
  #Calculate 7-day and 14-day rolling sum for daily cases and deaths
  mutate(cases_roll_7d  = rollsum(cases, 7, na.pad=TRUE, align="right")) %>% 
  mutate(death_roll_7d  = rollsum(death, 7, na.pad=TRUE, align="right")) %>% 
  mutate(cases_roll_14d = rollsum(cases, 14, na.pad=TRUE, align="right")) %>% 
  mutate(death_roll_14d = rollsum(death, 14, na.pad=TRUE, align="right")) %>% 
  #Generate 7-day lag variables for daily cases and deaths
  mutate(cases_roll_7d_lag7  = lag(cases_roll_7d,  n=7, default = NA)) %>% 
  mutate(cases_roll_14d_lag7 = lag(cases_roll_14d,  n=7, default = NA)) %>% 
  #Calculate CFR from the smoothed data with 0 and 7-day lags between the cases and deaths 
  mutate(CFR7_lag7      = death_roll_7d   / cases_roll_7d_lag7) %>%
  mutate(CFR14_lag7     = death_roll_14d  / cases_roll_14d_lag7) %>% 
  mutate(CFR7_same_day  = death_roll_7d   / cases_roll_7d) %>% 
  mutate(CFR14_same_day = death_roll_14d  / cases_roll_14d) %>% 
  ungroup() %>% 
  filter(week >= 20) #restrict to week 20 and after

### Identify outliers and missing CFR observations
# Should outliers be dropped (CFR > 1)?
outlier <- df2 %>% 
  filter(CFR7_lag7>1 | CFR14_lag7>1 | CFR7_same_day>1 | CFR14_same_day>1) %>% 
  select(countyname, state) %>% 
  distinct %>% 
  mutate(outlier.CFR = 1)

missing <- df2 %>% 
  filter(is.nan(CFR7_lag7) | is.nan(CFR14_lag7) | is.nan(CFR7_same_day) | is.nan(CFR14_same_day)) %>% 
  select(countyname, state) %>% 
  distinct %>% 
  mutate(miss.CFR = 1)


# If borth daily cases and deaths are 0 --> change CFR to 0 (R calculates it as NaN because the denominator is 0)
df2 <- df2 %>% 
  mutate(CFR7_lag7      = ifelse(death_roll_7d  == 0 & cases_roll_7d_lag7 == 0, 0, CFR7_lag7),
         CFR14_lag7     = ifelse(death_roll_14d == 0 & cases_roll_14d_lag7 == 0, 0, CFR14_lag7),
         CFR7_same_day  = ifelse(death_roll_7d  == 0 & cases_roll_7d == 0, 0, CFR7_same_day),
         CFR14_same_day = ifelse(death_roll_14d == 0 & cases_roll_14d == 0, 0, CFR14_same_day)) 
 
# Check again for missing obs.
missing <- df2 %>% 
  filter(is.nan(CFR7_lag7) | is.nan(CFR14_lag7) | is.nan(CFR7_same_day) | is.nan(CFR14_same_day)) %>% 
  select(countyname, state) %>% 
  distinct %>% 
  mutate(miss.CFR = 1)

## Note: Marin county has 1 Inf observation - because the denominator is 0 
inf <- df2 %>% 
  filter(CFR7_lag7==Inf | CFR14_lag7==Inf | CFR7_same_day==Inf | CFR14_same_day==Inf) %>% 
  select(countyname, state) %>% 
  distinct %>% 
  mutate(inf.CFR = 1)



### Conversion to daily (lara added)
df2_tot_covid <- df2 %>%                          
  group_by(countyname, state) %>%
  summarise(tot_covid = sum(death)) %>% 
  ungroup()


df2_daily <- df2 %>% 
  left_join(df2_tot_covid) %>% 
  filter(tot_covid >= 100) %>%                            #filter to counties with over than 100 deaths 
  #mutate(time = group_indices(.,day)) %>%                #NOTE: I changed "time" to start from 1 to last day in the period because one day was missing and creating a problem for the synth command in STATA          
  mutate(time = date-21914) %>%  
  mutate(heavy_smoke_p50  = ifelse(heavy_pct  >= 50, 1, 0)) %>% 
  mutate(heavy_smoke_p70  = ifelse(heavy_pct  >= 70, 1, 0)) %>%
  mutate(heavy_smoke_p90  = ifelse(heavy_pct  >= 90, 1, 0)) %>% 
  mutate(medium_smoke_p70 = ifelse(medium_pct >= 70, 1, 0)) %>% 
  mutate(id = group_indices(., countyname, state)) 


# Identify the start of the smoke
df2_daily_first_smoke_heavy50 <- df2_daily %>%           #find the first day with heavy smoke>50% 
  group_by(id, heavy_smoke_p50) %>% 
  filter(time == min(time)) %>% 
  slice(1) %>%
  ungroup() %>% 
  filter(heavy_smoke_p50==1) %>% 
  mutate(first_smoke_heavy50 = 1) %>% 
  select(countyname, state, time, first_smoke_heavy50)

df2_daily_first_smoke_heavy70 <- df2_daily %>%           #find the first day with heavy smoke>70% 
  group_by(id, heavy_smoke_p70) %>% 
  filter(time == min(time)) %>% 
  slice(1) %>%
  ungroup() %>% 
  filter(heavy_smoke_p70==1) %>% 
  mutate(first_smoke_heavy70 = 1) %>% 
  select(countyname, state, time, first_smoke_heavy70)

df2_daily_first_smoke_heavy90 <- df2_daily %>%           #find the first day with heavy smoke>90% 
  group_by(id, heavy_smoke_p90) %>% 
  filter(time == min(time)) %>% 
  slice(1) %>%
  ungroup() %>% 
  filter(heavy_smoke_p90==1) %>% 
  mutate(first_smoke_heavy90 = 1) %>% 
  select(countyname, state, time, first_smoke_heavy90)

df2_daily_first_smoke_medium70 <- df2_daily %>%           #find the first day with medium smoke>70% 
  group_by(id, medium_smoke_p70) %>% 
  filter(time == min(time)) %>% 
  slice(1) %>%
  ungroup() %>% 
  filter(medium_smoke_p70==1) %>% 
  mutate(first_smoke_medium70 = 1) %>% 
  select(countyname, state, time, first_smoke_medium70)


df2_daily <- df2_daily %>% 
  left_join(df2_daily_first_smoke_heavy50) %>% 
  mutate(smoke_heavy50_day1 = ifelse(first_smoke_heavy50 == 1, time, NA)) %>% 
  left_join(df2_daily_first_smoke_heavy70) %>% 
  mutate(smoke_heavy70_day1 = ifelse(first_smoke_heavy70 == 1, time, NA)) %>% 
  left_join(df2_daily_first_smoke_heavy90) %>% 
  mutate(smoke_heavy90_day1 = ifelse(first_smoke_heavy90 == 1, time, NA)) %>% 
  left_join(df2_daily_first_smoke_medium70) %>% 
  mutate(smoke_medium70_day1 = ifelse(first_smoke_medium70 == 1, time, NA)) %>% 
  filter(date >= 22067)     #filter data from beginning of June 


write.dta(df2_daily, "./Data/Full_data/Covid_smoke_by_county_daily_all_revised.dta", version = 12) # , label = attr(data, "label")



### Plot the smoothed data
fig2_1 <- df2_daily %>% filter(eligible == "Yes") %>% 
  filter(week >= 20) %>% 
  ggplot(aes(x=date2, y=death_roll_14d))+
  facet_wrap(~ countyname,  ncol = 6) + 
  theme(axis.title.x=element_blank())+
  geom_line()

fig2_2 <- df2_daily %>% filter(eligible == "Yes") %>% #Note: There is an unusual peak of cases for Marin County in July and August. Could this be due to misrecordsing?
  filter(week >= 20) %>% 
  ggplot(aes(x=date2, y=cases_roll_14d))+
  facet_wrap(~ countyname,  ncol = 6) + 
  theme(axis.title.x=element_blank())+
  geom_line()

fig2_3 <- df2_daily %>% filter(eligible == "Yes") %>%
  #filter(CFR14_lag7 < 0.15) %>%         #Note that I filter some extreme values. We should try investigate why these values are so high, possible misrecording 
  filter(week >= 20) %>%                      #restrict to week 20 and after
  ggplot(aes(x=date2, y=CFR14_lag7))+
  facet_wrap(~ countyname,  ncol = 6) + 
  theme(axis.title.x=element_blank())+
  geom_line()

fig2_4 <- df2_daily %>% filter(eligible == "Yes") %>%
  filter(CFR7_lag7 < 0.15) %>%              #Note that I filter some extreme values for Marin county; We should try investigate why these values are so high, possible misrecording 
  filter(week >= 20) %>%                        #restrict to week 20 and after
  ggplot(aes(x=date2, y=CFR7_lag7))+
  facet_wrap(~ countyname,  ncol = 6) + 
  theme(axis.title.x=element_blank())+
  geom_line()

fig2_5 <- df2_daily %>% filter(eligible == "Yes") %>%
  #filter(CFR7_same_day < 0.10) %>%            #Note that I filter some extreme values. We should try investigate why these values are so high, possible misrecording 
  filter(week >= 20) %>%         #restrict to week 20 to 40
  ggplot(aes(x=date2, y=CFR7_same_day))+
  facet_wrap(~ countyname,  ncol = 6) +
  theme(axis.title.x=element_blank())+
  geom_line()

fig2_6 <- df2_daily %>% filter(eligible == "Yes") %>%
  #filter(CFR14_same_day < 0.10) %>%            #Note that I filter some extreme values. We should try investigate why these values are so high, possible misrecording 
  filter(week >= 20) %>%         #restrict to week 20 to 40
  ggplot(aes(x=date2, y=CFR14_same_day))+
  facet_wrap(~ countyname,  ncol = 6) +
  theme(axis.title.x=element_blank())+
  geom_line()

fig2_7 <- grid.arrange(fig2_1, fig2_2, fig2_3, nrow=3, top = textGrob("Smooted data for eligible counties"))
ggsave("./Data/Diagnostics/Figures/Fig_2_Smoothed_data_14d.png", fig2_7, width = 22, height = 10, units = "in", dpi = 300)


fig2_8 <- grid.arrange(fig2_3, fig2_4, fig2_5, fig2_6, nrow=4, top = textGrob("Smooted data for eligible counties"))
ggsave("./Data/Diagnostics/Figures/Fig_3_Smoothed_data_alternative_CFR_measures.png", fig2_8, width = 22, height = 10, units = "in", dpi = 300)




# 3. Calculate weekly averages from the smoothed data
df2_weekly <- df2 %>% 
  group_by(countyname, state, week) %>%         
  summarise(CFR7_lag7     = mean(CFR7_lag7),         #calculate weekly averages 
            CFR14_lag7    = mean(CFR14_lag7),
            CFR7_same_day = mean(CFR7_same_day),
            CFR14_same_day= mean(CFR14_same_day), 
            heavy_pct     = mean(heavy_pct),
            medium_pct    = mean(medium_pct),
            death_roll_7d = mean(death_roll_7d),
            cases_roll_7d = mean(cases_roll_7d),
            mobility      = mean(populationstayingathome),
            mobility2     = mean(populationnotstayingathome)
            ) %>% 
  ungroup() %>% 
  left_join(df2_tot_covid) %>% 
  filter(tot_covid >= 100) %>%                            #filter to counties with over than 100 deaths 
  filter(week >= 20) %>%                                  #restrict to week 20 to 40
  mutate(time = week) %>%                           
  #mutate(smoke_bin = ifelse(heavy_pct >= 70, 1, 0)) %>% 
  mutate(heavy_smoke_p50  = ifelse(heavy_pct  >= 50, 1, 0)) %>% 
  mutate(heavy_smoke_p70  = ifelse(heavy_pct  >= 70, 1, 0)) %>%
  mutate(heavy_smoke_p90  = ifelse(heavy_pct  >= 90, 1, 0)) %>% 
  mutate(medium_smoke_p70 = ifelse(medium_pct >= 70, 1, 0)) %>% 
  mutate(id = group_indices(., countyname, state)) 


# Identify the start of the smoke
df2_first_smoke_heavy50 <- df2_weekly %>%           #find the first day with heavy smoke>50% 
  group_by(id, heavy_smoke_p50) %>% 
  filter(time == min(time)) %>% 
  slice(1) %>%
  ungroup() %>% 
  filter(heavy_smoke_p50==1) %>% 
  mutate(first_smoke_heavy50 = 1) %>% 
  select(countyname, state, time, first_smoke_heavy50)

df2_first_smoke_heavy70 <- df2_weekly %>%           #find the first day with heavy smoke>70% 
  group_by(id, heavy_smoke_p70) %>% 
  filter(time == min(time)) %>% 
  slice(1) %>%
  ungroup() %>% 
  filter(heavy_smoke_p70==1) %>% 
  mutate(first_smoke_heavy70 = 1) %>% 
  select(countyname, state, time, first_smoke_heavy70)

df2_first_smoke_heavy90 <- df2_weekly %>%           #find the first day with heavy smoke>90% 
  group_by(id, heavy_smoke_p90) %>% 
  filter(time == min(time)) %>% 
  slice(1) %>%
  ungroup() %>% 
  filter(heavy_smoke_p90==1) %>% 
  mutate(first_smoke_heavy90 = 1) %>% 
  select(countyname, state, time, first_smoke_heavy90)

df2_first_smoke_medium70 <- df2_weekly %>%           #find the first day with medium smoke>70% 
  group_by(id, medium_smoke_p70) %>% 
  filter(time == min(time)) %>% 
  slice(1) %>%
  ungroup() %>% 
  filter(medium_smoke_p70==1) %>% 
  mutate(first_smoke_medium70 = 1) %>% 
  select(countyname, state, time, first_smoke_medium70)


df2_weekly <- df2_weekly %>% 
  left_join(df2_first_smoke_heavy50) %>% 
  mutate(smoke_heavy50_week1 = ifelse(first_smoke_heavy50 == 1, time, NA)) %>% 
  left_join(df2_first_smoke_heavy70) %>% 
  mutate(smoke_heavy70_week1 = ifelse(first_smoke_heavy70 == 1, time, NA)) %>% 
  left_join(df2_first_smoke_heavy90) %>% 
  mutate(smoke_heavy90_week1 = ifelse(first_smoke_heavy90 == 1, time, NA)) %>% 
  left_join(df2_first_smoke_medium70) %>% 
  mutate(smoke_medium70_week1 = ifelse(first_smoke_medium70 == 1, time, NA)) %>% 
  mutate(eligible = ifelse( (countyname=="Alameda County"         |
                               countyname=="Contra Costa County"  | 
                               countyname=="Marin County"         |
                               countyname=="San Francisco County" |
                               countyname=="Santa Clara County"   |
                               countyname=="Sonoma County") &   state=="CA", "Yes", "No"))

write.dta(df2_weekly, "./Data/Full_data/Covid_smoke_by_county_week_all_revised.dta", version = 12) # , label = attr(data, "label")








