# Cox_Regression
유방암 여부와 근무 여부, 진통제 투여 정도와의 상관관계를 파악하기 위한 연구 논문 (분석 보조원으로 참여)

```r
### 유방암 환자 추출

brca <- c("C50","D05")

t20.c <- as.list(NA)
t20.c <- t20 %>% filter(substr(MAIN_SICK,1,3) %in% brca) %>% arrange(PERSON_ID, RECU_FR_DT)
length(unique(t20.c$PERSON_ID))     #4622

save(list=ls(),file="E:/seo_woman/age추가_fin/t20(~100).RData")

rm(list=ls())
load("E:/seo_woman/age추가_fin//t20(~100).RData")

# 정의 1. 1년 이내 3회 출현 - 외래(3876)
t20.cout <- t20.c %>% filter(FORM_CD == 3) %>% select(PERSON_ID, RECU_FR_DT)
t20.cout <- t20.cout %>% arrange(PERSON_ID, RECU_FR_DT)
length(unique(t20.cout$PERSON_ID))      #4570

brca_list <- split(t20.cout, t20.cout$PERSON_ID)

### 정의 1 함수 
# 15분 정도 걸림
case <- function(data,freq,term) {
  idlist <- c()
  n = freq-1 
  
  for(i in 1:length(data)){
    id <- data[[i]]
    list <- c()
    
    # freq 간격 내 기간이 term 이내인 환자를 return 함 
    #  ex. 3번의 암 진단을 받을 때까지의 기간이 365 일 이내이면 return -> 정의 1을 거꾸로 생각함
    for(j in 1:(nrow(id) - n)) {
      list[j] <- as.numeric(ymd(id$RECU_FR_DT[j+n]) - ymd(id$RECU_FR_DT[j]))
      if(sum(list <= term, na.rm = T) >= 1) {
        for(k in 1:length(list)) {
          idlist[[i]] <- id[k]
          if(list[k] <= term) break
        }
      }
    }
  }
  return(idlist)
}
# 1년(365) 이내에 암 코드가 3회 이상 출현한 환자  
bc1 <- case(brca_list, 3, 365)
bc1 <- bind_rows(bc1)

length(unique(bc1$PERSON_ID))

colnames(bc1) = c('PERSON_ID', 'CANCER_DATE')

# 정의 2. 입원 + 상병 (3364)
t20.cin <- t20.c %>% filter(FORM_CD == 2)
bc2 <- t20.cin %>% group_by(PERSON_ID) %>% summarise(CANCER_DATE = min(RECU_FR_DT))
length(unique(bc2$PERSON_ID))

# 정의 1, 2 합침 -> 두 가지 정의에 모두 해당하는 환자의 경우, 더 빠른 진단일자를 CANCER_DATE 로 (3987)
bc.final <- rbind(bc1, bc2) %>% 
  group_by(PERSON_ID) %>% 
  summarise(CANCER_DATE = min(CANCER_DATE))
```


## 2. 1:4 매칭
```r
### 일반인과 암환자 조인하여 매칭데이터 뽑기
matching_df<-rbind(df_can, df_cont)

if (!require(MatchIt)) install.packages('MatchIt')
library(MatchIt)

### 매칭하기 ################################
### exact는 오류가 발생하여 nearest로 함 ######
mat <- matchit(CANCER ~ RESIDENCE + IPSN_TYPE_CD, data = matching_df, method = "nearest", 
               exact = ~ RESIDENCE + IPSN_TYPE_CD, ratio = 4)
m0<-match.data(mat)
```

## cox-Regression
유방암 발생 여부와 발생까지의 기간을 종속변수로 고려하여 다양한 독립변수와의 상관관계를 파악
```r
table(df_d$CANCER)
table(df_d$work)

model1 <- coxph(Surv(duration, CANCER == 1) ~ work + AGE + income + 
                  BMI + alcohol + exer + 
                  CCI_CAT + nsaid, data = df_d)
summary(model1)
```

## ggsurvplot
누적 위험도에 관한 생존곡선 그래프를 일하는 여성과 일하지 않는 여성으로 나누어 작성
```r
# nonwork
ggsurvplot(survfit(Surv(duration/365, CANCER == 1) ~ nsaid, data=df_nw), data = df_nw, pval=FALSE, conf.int = TRUE, xlab = "Year",
           legend.labs = c("never", "non-regular", "regular"), legend.title = "Analgesics", linetype = "strata",
           risk.table = TRUE, risk.table.col = "strata", fun = "cumhaz", vpval.coord = c(0.1, 1), xlim = c(2,12), ylim = c(0,0.4), break.time.by = 2) 

# work
ggsurvplot(survfit(Surv(duration/365, CANCER == 1) ~ nsaid, data=df_w), data = df_w, pval=FALSE, conf.int = TRUE, xlab = "Year",
           legend.labs = c("never", "non-regular", "regular"), legend.title = "Analgesics", linetype = "strata",
           risk.table = TRUE, risk.table.col = "strata", fun = "cumhaz", vpval.coord = c(0.1, 1), xlim = c(2,12), ylim = c(0,0.4), break.time.by = 2) 

```
![image](https://user-images.githubusercontent.com/91238910/179354172-03db5feb-08bf-4fc5-a396-f7326ec0ee61.png)
