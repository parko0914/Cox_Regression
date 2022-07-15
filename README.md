# Cox_Regression
유방암 여부와 메트포민, 진통제 투여 정도와의 상관관계를 파악하기 위한 연구 논문 (분석 보조원으로 참여)

```r
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
