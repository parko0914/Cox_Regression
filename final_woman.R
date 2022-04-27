####Metformin Code
#options-------------------------------------------------------------
rm(list=ls())
options(scipen = 100) 
#R의 데이터를 올릴 수 있는 메모리 제한 해제
memory.limit((17592186044415))
gc() 
#---------------------------------------------------------------------

# Package -----------------------------------------------------------------
if (!require(dplyr)) install.packages('dplyr')
if (!require(comorbidity)) install.packages('comorbidity')
if (!require(data.table)) install.packages('data.table')
if (!require(reshape2)) install.packages('reshape2')
if (!require(lubridate)) install.packages('lubridate')
if (!require(bit64)) install.packages('bit64')
if (!require(readxl))  install.packages('readxl')
if (!require(finalfit)) install.packages('finalfit')
if (!require(ggplot2)) install.packages('ggplot2')
if (!require(survminer)) install.packages('survminer')
if (!require(survival)) install.packages('survival')
if (!require(Epi)) install.packages('Epi')
if (!require(broom)) install.packages('broom')
if (!require(plotrix)) install.packages('plotrix')
if (!require(gridExtra)) install.packages('gridExtra')
if (!require(ggjoy)) install.packages('ggjoy')
if (!require(naniar)) install.packages('naniar')
if (!require(psych)) install.packages('psych')
if (!require(epiDisplay)) install.packages('epiDisplay')
if (!require(oddsratio)) install.packages('oddsratio')
if (!require(GGally)) install.packages('GGally')
if (!require(ggfortify)) install.packages('ggfortify')
if (!require(magrittr)) install.packages('magrittr')
library(dplyr)
library(comorbidity)
library(data.table)
library(reshape2)
library(lubridate)
library(bit64)
library(readxl)
library(finalfit)
library(ggplot2)
library(survminer)
library(survival)
library(Epi)
library(broom)
library(plotrix)
library(gridExtra)
library(ggjoy)
library(naniar)
library(psych)
library(epiDisplay)
library(oddsratio)
library(GGally)
library(ggfortify)
library(magrittr)

select <- dplyr::select #dplyr 패키지 select함수 오류 방지

#############################################################
# 
# 진료DB: BC 정의 필요함                                 
#  - BC 정의: 1년 내 C50 / D05 출현 3회 또는 C/D + 입원
# 자격DB: BC와 Non-BC 자료 필요함
#  - 매칭 변수: age / area 
# 매칭: Exact 또는 Nearest
# 
#############################################################

# 매칭 데이터셋 구축 순서
# 1. BC 환자 추출 - BC 자격 / 검진 데이터셋 
# (1) BC 환자 조작적 정의 : 1년 내 암 코드 3회 출현 또는 입원
# (2) BC 환자 추출 뒤 암 진단 연도와 ID 리스트 
# (3) BC 환자 자격 / 검진 자료 
# 2. Non-BC 추출 - 자격 / 검진 데이터셋
# (1) BC 환자 자격 / 검진 자료 
# (2) Non-BC 환자 자격 / 검진 자료 
# 3. 매칭 데이터셋 구축
# (1) 구축한 자료 합치기
# 4. 매칭

# 데이터 불러오기 (t20) : 유방암 환자 추출

setwd('E:\\edata\\medicaldata\\20t')
listfile <- list.files() # 폴더 내 파일명 저장
t20 <- list() # 빈 리스트 생성 뒤 연도별 데이터 input

# 필요한 변수만 추출
# 진료DB 명세서(t20)에서 
# 개인일련번호(PERSON_ID), 
# 청구일련번호(KEY_SEQ), 
# 요양개시일자(RECU_FR_DT), 
# FORM_CD,
# 주상병(MAIN_SICK) 및 부상병(SUB_SICK)
# 총처방일수(TOT_PRES_DD_CNT)


for(i in 1:12){t20[[i]] <- data.table::fread(listfile[i],
                                             stringsAsFactors = F,
                                             encoding = "Latin-1") %>% 
  dplyr::select(PERSON_ID, 
                KEY_SEQ, 
                RECU_FR_DT, 
                FORM_CD,
                MAIN_SICK,
                SUB_SICK,
                TOT_PRES_DD_CNT)
}

t20 <- bind_rows(t20)

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

save(t20.c,bc.final, bc1,bc2 ,file="E:/seo_woman/age추가_fin/bc_new.RData")

#############################################################################
### baseline 찾기
rm(list=ls())
load("E:/seo_woman/age추가_fin/bc_new.RData")

#######################
### 자격 DB ###
########################
path <- "E:/edata/medicaldata/jk/"
filelist <- list.files(path)
jk <- list()
for (i in 1:length(filelist)) {
  jk[[i]] <- fread(paste0(path, filelist[i]), 
                   stringsAsFactors = F, 
                   encoding = "Latin-1") %>% 
    dplyr::select(PERSON_ID, 
                  STND_Y, 
                  SEX, AGE_GROUP, SIDO,
                  CTRB_PT_TYPE_CD,
                  IPSN_TYPE_CD,
                  DTH_YM)
}
rm(filelist, i, path)

jk <- dplyr::bind_rows(jk)

jk <- jk %>% arrange(PERSON_ID, STND_Y) %>% 
  filter(PERSON_ID %in% bc.final$PERSON_ID)

### 첫 진입 날짜 기준 jk 데이터 df 에 결합
jk_first <- jk %>% group_by(PERSON_ID) %>% summarise(firstyear = nth(STND_Y, 1))
bc.final <- bc.final %>% left_join(jk_first)

jk <- jk %>% rename(firstyear = STND_Y)
bc.final <- bc.final %>% left_join(jk)

##############################
### 검진 DB 추가 ###
##############################
path <- "E:\\edata\\medicaldata\\gj\\"
filelist <- list.files(path)
gj <- list()
for (i in 1:7) {
  gj[[i]] <- fread(paste0(path, filelist[i]), 
                   stringsAsFactors = F, 
                   encoding = "Latin-1") %>% 
    mutate(WAIST = NA,
           TRIGLYCERIDE = NA,
           HDL_CHOLE = NA,
           LDL_CHOLE = NA,
           CREATININE = NA,
           MOV20_WEK_FREQ_ID = NA,
           MOV30_WEK_FREQ_ID = NA,
           WLK30_WEK_FREQ_ID = NA,
           HCHK_HPRTS_PMH_YN = NA,
           HCHK_DIABML_PMH_YN = NA) %>% 
    dplyr::select(PERSON_ID, 
                  HCHK_YEAR, 
                  HEIGHT, WEIGHT, 
                  WAIST, 
                  BP_HIGH, BP_LWST,
                  BLDS, 
                  TOT_CHOLE, TRIGLYCERIDE, HDL_CHOLE, LDL_CHOLE, 
                  HMG, 
                  OLIG_PROTE_CD, 
                  CREATININE, SGOT_AST, SGPT_ALT, GAMMA_GTP, 
                  SMK_STAT_TYPE_RSPS_CD, DRNK_HABIT_RSPS_CD,
                  EXERCI_FREQ_RSPS_CD,
                  HCHK_PMH_CD1,HCHK_PMH_CD2,HCHK_PMH_CD3)
}
for (i in 8:length(filelist)) {
  gj[[i]] <- fread(paste0(path, filelist[i]), 
                   stringsAsFactors = F, 
                   encoding = "Latin-1") %>% 
    mutate(EXERCI_FREQ_RSPS_CD = NA,
           HCHK_PMH_CD1 = NA,
           HCHK_PMH_CD2 = NA,
           HCHK_PMH_CD1 = NA) %>% 
    dplyr::select(PERSON_ID, 
                  HCHK_YEAR, 
                  HEIGHT, WEIGHT, 
                  WAIST, 
                  BP_HIGH, BP_LWST,
                  BLDS, 
                  TOT_CHOLE, TRIGLYCERIDE, HDL_CHOLE, LDL_CHOLE, 
                  HMG, 
                  OLIG_PROTE_CD, 
                  CREATININE, SGOT_AST, SGPT_ALT, GAMMA_GTP, 
                  SMK_STAT_TYPE_RSPS_CD, DRNK_HABIT_RSPS_CD,
                  MOV20_WEK_FREQ_ID, MOV30_WEK_FREQ_ID, WLK30_WEK_FREQ_ID, 
                  HCHK_HPRTS_PMH_YN, HCHK_DIABML_PMH_YN)
}
rm(filelist, i, path)

gj <- bind_rows(gj) %>% 
  arrange(PERSON_ID, HCHK_YEAR) %>% 
  filter(PERSON_ID %in% bc.final$PERSON_ID)
length(unique(gj$PERSON_ID))                  #2923

gj.final <- gj %>% group_by(PERSON_ID) %>% summarise(FIRST_HCHK = nth(HCHK_YEAR,1))
bc.final <- bc.final %>% left_join(gj.final, by='PERSON_ID')

# 제외 기준 정리
### 여성만 추출 
# 1. 남성 제외 : 18명, 제외 후 3969명
bc.man <- bc.final %>% filter(SEX==1)
bc.final <- bc.final %>% filter(SEX==2)

# 2. 첫검진이 암발생일 뒤에 있는 경우 제외 : 814명 제외 후, 3155 명
bc.first_c <- bc.final %>% dplyr::filter(FIRST_HCHK > year(lubridate::ymd(CANCER_DATE)))
bc.final <- bc.final %>% dplyr::filter(!PERSON_ID %in% bc.first_c$PERSON_ID)

### 연령대 변수 추가 
# 3-1. 나이 제외 기준 (25세 미만 환자 제외) : 60 명
under25 <- bc.final %>% filter(AGE_GROUP < 6)
# 3-2. 나이 제외 기준 (65세 이상 환자 제외) : 228 명
up65 <- bc.final %>% filter(AGE_GROUP > 13)
# 최종 df 나이 그룹 : 1. 25-39 (6,7,8) 2. 40-49 (9,10) 3. 50-64 (11,12,13) : 2867 명 
bc.final<- bc.final %>% filter(AGE_GROUP %in% c(6,7,8,9,10,11,12,13))

### 4. 첫 검진년도 2002년인 환제 제외 : 449명 제외 후, 2418명
gj2002 <- bc.final %>% dplyr::filter(FIRST_HCHK == 2002)
bc.final <- bc.final %>% filter(!PERSON_ID %in% gj2002$PERSON_ID)

### 5. 의료 program 대상자 제외 : 50명, 이 대상자들은 추출 후 제외함
medi_re <- bc.final %>% filter(IPSN_TYPE_CD %in% c(7,8))

### 6. 검진 자료 없는 환자 제외 : 910명 제외 후, 1058 명
gjno <- bc.final %>% filter(is.na(FIRST_HCHK))

df.final <- bc.final
save(list=ls() ,file="E:/seo_woman/age추가_fin/bc_df_final.RData")


### 암환자 데이터에서 필요한 변수들만 뽑아서 정리##############
### PERSON_ID, firstyear, 연령대, 거주지, CANCER #############
###############################################################
rm(list=ls())
load("E:/seo_woman/age추가_fin/bc_df_final.RData")

matching_df<- df.final %>% dplyr::select(PERSON_ID, AGE_GROUP, SIDO, IPSN_TYPE_CD, firstyear)

### 거주지 변수 추가 (수도권, 비수도권)
### 서울, 인천, 경기도면 0, 나머지는 1
matching_df<- matching_df %>% mutate(RESIDENCE = ifelse(SIDO == 11 | SIDO == 41 | SIDO == 28, 0, 1))

### 기존 SIDO 변수 제거
matching_df$SIDO <- NULL

### CANCER 변수 추가
matching_df$CANCER <- 1

### 매칭을 위한 암환자 df
df_can<-matching_df

#########################################################################
### 암환자가 아닌 일반인 환자 뽑기 ######################################
#########################################################################
### 우선 자격 데이터 에서 PERSON_ID와 STND_Y 을 뽑기 (baseline 을 코호트 첫 진입으로 변경)
path <- "E:/edata/medicaldata/jk/"
filelist <- list.files(path)
jk <- list()
for (i in 1:length(filelist)) {
  jk[[i]] <- fread(paste0(path, filelist[i]), 
                   stringsAsFactors = F, 
                   encoding = "Latin-1") %>% 
    dplyr::select(PERSON_ID, 
                  STND_Y, 
                  SEX, AGE_GROUP, SIDO,
                  CTRB_PT_TYPE_CD,
                  IPSN_TYPE_CD,
                  DTH_YM)
}
rm(filelist, i, path)

jk <- dplyr::bind_rows(jk)

jk_m <- jk %>% arrange(PERSON_ID, STND_Y) %>% 
  filter(!PERSON_ID %in% df_can$PERSON_ID)
length(unique(jk_m$PERSON_ID))          #1124183

### 첫 진입년도 뽑기 ###################################################
df_cont <- jk_m %>% group_by(PERSON_ID) %>% summarise(firstyear = nth(STND_Y,1))
rm(jk_m)

### 자격 데이터 추가하기 전에 검진 자료 없는 환자는 미리 빼자
path <- "E:\\edata\\medicaldata\\gj\\"
filelist <- list.files(path)
gj <- list()
for (i in 1:7) {
  gj[[i]] <- fread(paste0(path, filelist[i]), 
                   stringsAsFactors = F, 
                   encoding = "Latin-1") %>% 
    mutate(WAIST = NA,
           TRIGLYCERIDE = NA,
           HDL_CHOLE = NA,
           LDL_CHOLE = NA,
           CREATININE = NA,
           MOV20_WEK_FREQ_ID = NA,
           MOV30_WEK_FREQ_ID = NA,
           WLK30_WEK_FREQ_ID = NA,
           HCHK_HPRTS_PMH_YN = NA,
           HCHK_DIABML_PMH_YN = NA) %>% 
    dplyr::select(PERSON_ID, 
                  HCHK_YEAR, 
                  HEIGHT, WEIGHT, 
                  WAIST, 
                  BP_HIGH, BP_LWST,
                  BLDS, 
                  TOT_CHOLE, TRIGLYCERIDE, HDL_CHOLE, LDL_CHOLE, 
                  HMG, 
                  OLIG_PROTE_CD, 
                  CREATININE, SGOT_AST, SGPT_ALT, GAMMA_GTP, 
                  SMK_STAT_TYPE_RSPS_CD, DRNK_HABIT_RSPS_CD,
                  EXERCI_FREQ_RSPS_CD,
                  HCHK_PMH_CD1,HCHK_PMH_CD2,HCHK_PMH_CD3)
}
for (i in 8:length(filelist)) {
  gj[[i]] <- fread(paste0(path, filelist[i]), 
                   stringsAsFactors = F, 
                   encoding = "Latin-1") %>% 
    mutate(EXERCI_FREQ_RSPS_CD = NA,
           HCHK_PMH_CD1 = NA,
           HCHK_PMH_CD2 = NA,
           HCHK_PMH_CD1 = NA) %>% 
    dplyr::select(PERSON_ID, 
                  HCHK_YEAR, 
                  HEIGHT, WEIGHT, 
                  WAIST, 
                  BP_HIGH, BP_LWST,
                  BLDS, 
                  TOT_CHOLE, TRIGLYCERIDE, HDL_CHOLE, LDL_CHOLE, 
                  HMG, 
                  OLIG_PROTE_CD, 
                  CREATININE, SGOT_AST, SGPT_ALT, GAMMA_GTP, 
                  SMK_STAT_TYPE_RSPS_CD, DRNK_HABIT_RSPS_CD,
                  MOV20_WEK_FREQ_ID, MOV30_WEK_FREQ_ID, WLK30_WEK_FREQ_ID, 
                  HCHK_HPRTS_PMH_YN, HCHK_DIABML_PMH_YN)
}
rm(filelist, i, path)

gj <- bind_rows(gj)
length(unique(gj$PERSON_ID))

# 첫 검진 년도가 2002년인 환자 뽑아서 제거해주고,
gj_2002 <- gj %>% group_by(PERSON_ID) %>% summarise(fy = min(HCHK_YEAR))
gj_2002 <- gj_2002 %>% filter(fy == 2002)

gj <- gj %>% filter(!PERSON_ID %in% gj_2002$PERSON_ID)

# gj 자료(!= 2002)가 있는 df_cont
df_cont <- df_cont %>% filter(PERSON_ID %in% gj$PERSON_ID)

### df_cont에 자격 데이터 추가 ###########################################
### 1.1.3 자격 DB
path <- "E:/edata/medicaldata/jk/"
filelist <- list.files(path)
jk <- list()
for (i in 1:length(filelist)) {
  jk[[i]] <- fread(paste0(path, filelist[i]), 
                   stringsAsFactors = F, 
                   encoding = "Latin-1") %>% 
    dplyr::select(PERSON_ID, 
                  STND_Y, 
                  SEX, AGE_GROUP, SIDO,
                  CTRB_PT_TYPE_CD,
                  IPSN_TYPE_CD,
                  DTH_YM)
}
rm(filelist, i, path)

jk <- dplyr::bind_rows(jk)

jk <- jk %>% arrange(PERSON_ID, STND_Y) %>% 
  filter(PERSON_ID %in% df_cont$PERSON_ID)

### 여자, 25~64세만 뽑기 ###############################################
### 처음 1124385 명
jk <- jk %>% filter(SEX==2)           
jk <- subset(jk, AGE_GROUP >= 6)   
jk <- subset(jk, AGE_GROUP < 14)    

### 변수 정리 ##########################################################
### 1. jk와 df_cont을 join할때 PERSON_ID와 FIRST_HCHK을 해줘야하므로 이름 변경
jk$firstyear<-jk$STND_Y

### 거주지 변수 추가
### 서울, 인천, 경기도면 0, 나머지는 1
jk<- jk %>% mutate(RESIDENCE = ifelse(SIDO == 11 | SIDO == 41 | SIDO == 28, 0, 1))

### CANCER변수 추가
jk$CANCER <- 0

### 필요없는 변수 제거 #################################################
jk$STND_Y<-NULL
jk$SEX<-NULL
jk$SIDO<-NULL
jk$CTRB_PT_TYPE_CD<-NULL
jk$DTH_YM<-NULL

### 조인 하기전 df_cont을 jk에 있는 PERSON_ID만 남아 있게 필터 
df_cont <- df_cont %>% filter(PERSON_ID  %in% jk$PERSON_ID)
length(unique(jk$PERSON_ID))        #228789
length(unique(df_cont$PERSON_ID))   #228789

### 조인하기
df_cont <- df_cont %>% left_join(jk)

### NA 제거
df_cont <-na.omit(df_cont)

### 최종 일반인 : 185863 명
### 최종 암환자 : 1508명
save(list=ls() ,file="E:/seo_woman/age추가_fin/bc_control_final.RData")

#########################################################################
############### 매칭 데이터셋 뽑기 ############################
rm(list=ls())
load("E:/seo_woman/age추가_fin/bc_control_final.RData")

### 일반인과 암환자 조인하여 매칭데이터 뽑기
matching_df<-rbind(df_can, df_cont)

if (!require(MatchIt)) install.packages('MatchIt')
library(MatchIt)

### 매칭하기 ################################
### exact는 오류가 발생하여 nearest로 함 ######
mat <- matchit(CANCER ~ firstyear + AGE_GROUP + RESIDENCE + IPSN_TYPE_CD, method = "nearest", data = matching_df, ratio = 4)
m<-match.data(mat)

# 확인 작업 
chisq.test(m$AGE_GROUP, m$CANCER)
chisq.test(m$RESIDENCE, m$CANCER)
chisq.test(m$IPSN_TYPE_CD, m$CANCER)

length(which(m$CANCER == 0))

xtabs(~AGE_GROUP + CANCER, data = m)

save(m, file="E:/seo_woman/age추가_fin/final_df_new.RData")

################################################################################
###################### 여기까지가 매칭 df 뽑는 과정 ############################
## 아래 부터는 1. bc_df_final.RData (암환자만 있는 데이터) ####################
## 2. final_df_new.RData (case, control 군 둘다 포함) 만 사용하면 됨 ##########
############## (재매칭하게 되면 데이터셋 다시 바뀌므로) ########################
################################################################################

##################################################################
### 최종 데이터 셋 (m 데이터) ##############################
### 이후 이 데이터셋에 해당하는 PERSON_ID를 이용하여 
### 20t, gj, jk는 PERSON_ID로
### 30t, 40t, 60t는 20t의 KEY_SEQ를 이용하여 필터
rm(list=ls())
load("E:/seo_woman/age추가_fin/final_df_new.RData")
load("E:/seo_woman/age추가_fin/bc_df_final.RData")

select <- dplyr::select

#### 기존에 있던 암 환자 데이터에서 CANCER_DATE 만 불러와서 붙이기
df <- m
rm(m)

df <- df %>% left_join(df.final %>% dplyr::select(PERSON_ID, CANCER_DATE))
rm(df.final, gj,jk,jk_first,t20.c, bc.final, bc2)

# AGE 생성
## 25-39 0, 40대 1, 50-64 2
df<- df %>% mutate(AGE = ifelse(AGE_GROUP>=6 & AGE_GROUP<9, 0,
                                ifelse(AGE_GROUP>=9 & AGE_GROUP<11,1,
                                       ifelse(AGE_GROUP>=11 & AGE_GROUP<14,2,NA))))

### jk 불러와서 insurance(자격 보험), income(소득) 생성
setwd('E:/edata/medicaldata/jk')
listfile <- list.files()
jk <- list()
for(i in 1:12){jk[[i]] <- data.table::fread(listfile[i]) %>% 
  dplyr::filter(PERSON_ID %in% df$PERSON_ID)}
jk <- dplyr::bind_rows(jk)
length(unique(jk$PERSON_ID))

### firstyear 맞춰주고 붙이기
jk$firstyear <- jk$STND_Y
df <- df %>% left_join(jk %>% select(c(PERSON_ID, firstyear, CTRB_PT_TYPE_CD)))

df <- df %>% mutate(insurance = ifelse(IPSN_TYPE_CD %in% c(1,2), 0,
                                       ifelse(IPSN_TYPE_CD %in% c(5,6), 1, 2)),
                    income = ifelse(CTRB_PT_TYPE_CD == 0, 0,
                                    ifelse(CTRB_PT_TYPE_CD %in% c(1,2,3), 1,
                                           ifelse(CTRB_PT_TYPE_CD %in% c(4,5,6), 2, 3))))

xtabs(~df$income + df$CANCER)
xtabs(~df$insurance + df$CANCER)

### 검진 자료 관련 변수 붙이기
path <- "E:\\edata\\medicaldata\\gj\\"
filelist <- list.files(path)
gj <- list()
for (i in 1:7) {
  gj[[i]] <- fread(paste0(path, filelist[i]), 
                   stringsAsFactors = F, 
                   encoding = "Latin-1") %>% 
    mutate(WAIST = NA,
           TRIGLYCERIDE = NA,
           HDL_CHOLE = NA,
           LDL_CHOLE = NA,
           CREATININE = NA,
           MOV20_WEK_FREQ_ID = NA,
           MOV30_WEK_FREQ_ID = NA,
           WLK30_WEK_FREQ_ID = NA,
           HCHK_HPRTS_PMH_YN = NA,
           HCHK_DIABML_PMH_YN = NA) %>% 
    dplyr::select(PERSON_ID, 
                  HCHK_YEAR, 
                  HEIGHT, WEIGHT, 
                  WAIST, 
                  BP_HIGH, BP_LWST,
                  BLDS, 
                  TOT_CHOLE, TRIGLYCERIDE, HDL_CHOLE, LDL_CHOLE, 
                  HMG, 
                  OLIG_PROTE_CD, 
                  CREATININE, SGOT_AST, SGPT_ALT, GAMMA_GTP, 
                  SMK_STAT_TYPE_RSPS_CD, DRNK_HABIT_RSPS_CD,
                  EXERCI_FREQ_RSPS_CD,
                  HCHK_PMH_CD1,HCHK_PMH_CD2,HCHK_PMH_CD3)
}
for (i in 8:length(filelist)) {
  gj[[i]] <- fread(paste0(path, filelist[i]), 
                   stringsAsFactors = F, 
                   encoding = "Latin-1") %>% 
    mutate(EXERCI_FREQ_RSPS_CD = NA,
           HCHK_PMH_CD1 = NA,
           HCHK_PMH_CD2 = NA,
           HCHK_PMH_CD1 = NA) %>% 
    dplyr::select(PERSON_ID, 
                  HCHK_YEAR, 
                  HEIGHT, WEIGHT, 
                  WAIST, 
                  BP_HIGH, BP_LWST,
                  BLDS, 
                  TOT_CHOLE, TRIGLYCERIDE, HDL_CHOLE, LDL_CHOLE, 
                  HMG, 
                  OLIG_PROTE_CD, 
                  CREATININE, SGOT_AST, SGPT_ALT, GAMMA_GTP, 
                  SMK_STAT_TYPE_RSPS_CD, DRNK_HABIT_RSPS_CD,
                  MOV20_WEK_FREQ_ID, MOV30_WEK_FREQ_ID, WLK30_WEK_FREQ_ID, 
                  HCHK_HPRTS_PMH_YN, HCHK_DIABML_PMH_YN)
}
rm(filelist, i, path)


gj <- bind_rows(gj) %>% filter(PERSON_ID %in% df$PERSON_ID)

### 첫 검진날짜에 해당하는 자료 뽑아와서 붙이기
gj1 <- gj %>% group_by(PERSON_ID) %>% summarise(FIRST_YEAR = min(HCHK_YEAR))
gj <- gj %>% left_join(gj1) %>% filter(HCHK_YEAR == FIRST_YEAR)
rm(gj1)

# BMI
gj <- gj %>% mutate(BMI = WEIGHT / (HEIGHT/100)**2)

# Alcohol
gj <- gj %>% 
  mutate(Alcohol = ifelse(HCHK_YEAR >= 2009 & DRNK_HABIT_RSPS_CD ==1, 0,
                          ifelse(HCHK_YEAR >= 2009 & DRNK_HABIT_RSPS_CD <8, 1,
                                 ifelse(HCHK_YEAR >= 2009 & DRNK_HABIT_RSPS_CD ==8, 2,
                                        ifelse(HCHK_YEAR < 2009 & DRNK_HABIT_RSPS_CD ==1, 0,
                                               ifelse(HCHK_YEAR < 2009 & DRNK_HABIT_RSPS_CD <5, 1,
                                                      ifelse(HCHK_YEAR < 2009 & DRNK_HABIT_RSPS_CD ==5, 2, NA)))))))
# Smoking
gj$Smoking <- gj$SMK_STAT_TYPE_RSPS_CD
table(gj$Smoking)

# Exercise
gj <- gj %>% 
  mutate(Exercise = ifelse(HCHK_YEAR >= 2009 & (MOV20_WEK_FREQ_ID %in% c(7,8) | MOV30_WEK_FREQ_ID %in% c(7,8) | WLK30_WEK_FREQ_ID %in% c(7,8)), 2,
                           ifelse(HCHK_YEAR >= 2009 & (MOV20_WEK_FREQ_ID %in% c(2,3,4,5,6) | MOV30_WEK_FREQ_ID %in% c(2,3,4,5,6) | WLK30_WEK_FREQ_ID %in% c(2,3,4,5,6)), 1,
                                  ifelse(HCHK_YEAR >= 2009 & MOV20_WEK_FREQ_ID == 1 & MOV30_WEK_FREQ_ID == 1 & WLK30_WEK_FREQ_ID == 1, 0,
                                         ifelse(HCHK_YEAR < 2009 & EXERCI_FREQ_RSPS_CD == 5, 2,
                                                ifelse(HCHK_YEAR < 2009 & EXERCI_FREQ_RSPS_CD %in% c(2,3,4), 1,
                                                       ifelse(HCHK_YEAR < 2009 & EXERCI_FREQ_RSPS_CD == 1, 0, NA)))))))
table(gj$Exercise)

df <- df %>% left_join(gj %>% select(PERSON_ID, BMI, Alcohol, Smoking, Exercise, TOT_CHOLE))

### 검진변수 범주화
table(df$Alcohol)
table(df$Smoking)
table(df$Exercise)
sum(is.na(df$Alcohol))
sum(is.na(df$Smoking))
sum(is.na(df$Exercise))

df <- df %>% mutate(bmi = ifelse(BMI < 25, 0, 1),
                    alcohol = ifelse(Alcohol == 0,0,1),
                    smoke = ifelse(Smoking == 1, 0, 1),
                    exer = ifelse(Exercise == 0, 0, 1))

### 고혈압, 당뇨 정의를 위해 t20, t60 불러오기
## t20
setwd('E:\\edata\\medicaldata\\20t')
listfile <- list.files()
t20 <- list()
for(i in 1:12){t20[[i]] <- data.table::fread(listfile[i],
                                             stringsAsFactors = F,
                                             encoding = "Latin-1") %>% 
  dplyr::select(PERSON_ID, 
                KEY_SEQ, 
                RECU_FR_DT, 
                MAIN_SICK, 
                SUB_SICK, 
                TOT_PRES_DD_CNT,
                FORM_CD) %>% 
  dplyr::filter(PERSON_ID %in% df$PERSON_ID)}
t20 <- dplyr::bind_rows(t20)


## t60 
setwd('E:\\edata\\medicaldata\\60t')
listfile <- list.files()
t60 <- list()
for(i in 1:12){t60[[i]] <- data.table::fread(listfile[i]) %>% 
  dplyr::select(KEY_SEQ, 
                RECU_FR_DT,
                GNL_NM_CD,
                MDCN_EXEC_FREQ) %>% 
  dplyr::filter(KEY_SEQ %in% t20$KEY_SEQ)}
t60 <- dplyr::bind_rows(t60)


save(list=ls(), file="E:/seo_woman/age추가_fin/newbaseline_jk_gj_t2060.RData")

rm(list=ls())
load("E:/seo_woman/age추가_fin/newbaseline_jk_gj_t2060.RData")

########## 당뇨환자와 고혈압 환자 정의하기 ##########
### hypertension 기준 
# 기준 1 : 첫 검진 년도 전까지의 진료기록 중에서 주/부상병에 고혈압 코드 2회이상 발견
t20_hyper <- t20 %>% inner_join(gj %>% select(PERSON_ID, HCHK_YEAR)) %>% 
  filter(year(ymd(RECU_FR_DT)) <= HCHK_YEAR) %>% arrange(PERSON_ID, RECU_FR_DT)
t20_hyper <- t20_hyper %>% filter(substr(MAIN_SICK,1,3) %in% c('I10', 'I11', 'I12', 'I13', 'I15') |
                                    substr(SUB_SICK,1,3) %in% c('I10', 'I11', 'I12', 'I13', 'I15'))
hyper_1 <- t20_hyper %>% group_by(PERSON_ID) %>% summarise(count= n())
hyper_1 <- hyper_1 %>% filter(count>=2)
length(unique(hyper_1$PERSON_ID))


# 기준 2 : 고혈압 코드와 당뇨 약제 동시에 존재
medicine_h <- fread('E:/seo_woman/0_인수인계 자료입니다!_210619/hypertension.csv', 
                    stringsAsFactors = F, 
                    encoding = "Latin-1")
medicine_h <- medicine_h %>% 
  dplyr::mutate(ATC_CAT = substr(code, 1, 4)) %>% 
  unique()
t60_hyper <- t60 %>% 
  dplyr::mutate(ATC_CAT = substr(GNL_NM_CD, 1, 4)) %>% 
  dplyr::filter(ATC_CAT %in% medicine_h$ATC_CAT) %>% 
  dplyr::inner_join(medicine_h) %>% 
  dplyr::left_join(t20 %>% 
                     dplyr::select(PERSON_ID,
                                   KEY_SEQ))

# 진료기록에서 고혈압 코드가 있는 환자들 필터
hyper_2 <- t60_hyper %>% filter(PERSON_ID %in% t20_hyper$PERSON_ID)
length(unique(hyper_2$PERSON_ID))

hyper <- c(hyper_1$PERSON_ID, hyper_2$PERSON_ID) %>% unique()

df <- df %>% mutate(hypertension = ifelse(PERSON_ID %in% hyper, 1, 0))



### diabetes 기준
# 기준 1 : 첫 검진 년도 전까지의 진료기록 중에서 주/부상병에 당뇨 코드 2회이상 발견
t20_d <- t20 %>% inner_join(gj %>% select(PERSON_ID, HCHK_YEAR)) %>% 
  filter(year(ymd(RECU_FR_DT)) <= HCHK_YEAR) %>% arrange(PERSON_ID, RECU_FR_DT)
t20_d <- t20_d %>% filter(substr(MAIN_SICK,1,3) %in% c('E11', 'E12', 'E13', 'E14') |
                            substr(SUB_SICK,1,3) %in% c('E11', 'E12', 'E13', 'E14'))
d_1 <- t20_d %>% group_by(PERSON_ID) %>% summarise(count= n())
d_1 <- d_1 %>% filter(count>=2)
length(unique(d_1$PERSON_ID))

# 기준 2 : 당뇨 코드와 당뇨 약제 동시에 존재
medicine_d <- fread('E:/seo_woman/0_인수인계 자료입니다!_210619/당뇨약제.csv', 
                    stringsAsFactors = F, 
                    encoding = "Latin-1")
medicine_d <- medicine_d %>% 
  dplyr::mutate(ATC_CAT = substr(code, 1, 4)) %>% 
  unique()

t60_d <- t60 %>% 
  dplyr::mutate(ATC_CAT = substr(GNL_NM_CD, 1, 4)) %>% 
  dplyr::filter(ATC_CAT %in% medicine_d$ATC_CAT) %>% 
  dplyr::inner_join(medicine_d) %>% 
  dplyr::left_join(t20 %>% 
                     dplyr::select(PERSON_ID,
                                   KEY_SEQ))

# 진료기록에서 당뇨 코드가 있는 환자들 필터
d_2 <- t60_d %>% filter(PERSON_ID %in% t20_d$PERSON_ID)
length(unique(d_2$PERSON_ID))

d_t <- c(d_1$PERSON_ID, d_2$PERSON_ID) %>% unique()
df <- df %>% mutate(diabetes = ifelse(PERSON_ID %in% d_t, 1, 0))

rm(hyper_1, hyper_2, d_1, d_2, t20_hyper, t60_hyper, t20_d, t60_d)

### childbirth (출산 여부) 기준
# t20 에서 주/부상병 첫 글자가 O 이면 ever, 아니면 never
# 단, 암 환자의 경우, 암이 발생하기 전까지 기간만 고려함
t20_c <- t20 %>% filter(substr(MAIN_SICK,1,1) == "O" | 
                          substr(SUB_SICK,1,1) == "O")
t20_c <- t20_c %>% inner_join(df %>% select(PERSON_ID, CANCER_DATE))

t20_c1 <- t20_c %>% filter(is.na(CANCER_DATE) == F) %>% 
  filter(ymd(RECU_FR_DT) <= ymd(CANCER_DATE))
t20_c0 <- t20_c %>% filter(is.na(CANCER_DATE))
t20_c <- rbind(t20_c0, t20_c1) %>% unique()

df <- df %>% mutate(childbirth = ifelse(PERSON_ID %in% t20_c$PERSON_ID,1,0))
table(df$childbirth)

### CCI 변수 추가 
t20 <- t20 %>% filter(PERSON_ID %in% df$PERSON_ID)

setwd('E:\\edata\\medicaldata\\40t')
listfile <- list.files()
t40 <- list()
for(i in 1:12){t40[[i]] <- data.table::fread(listfile[i]) %>% 
  dplyr::select(KEY_SEQ,
                RECU_FR_DT,
                SICK_SYM) %>% 
  dplyr::filter(KEY_SEQ %in% t20$KEY_SEQ)}

t40 <- dplyr::bind_rows(t40) %>% 
  dplyr::inner_join(t20 %>% 
                      dplyr::select(PERSON_ID,
                                    KEY_SEQ)) %>% 
  dplyr::filter(PERSON_ID %in% df$PERSON_ID) 

# t40에 firstyear 붙여서, 코호트 진입 이후 1년까지의 데이터 필터
t40_CCI <- t40 %>% left_join(df %>% select(PERSON_ID, firstyear))
t40_CCI <- t40_CCI %>% filter(year(ymd(RECU_FR_DT)) <= firstyear + 1)
length(unique(t40_CCI$PERSON_ID))

score <- comorbidity::comorbidity(t40_CCI,
                                  id = 'PERSON_ID',
                                  code = 'SICK_SYM',
                                  score = 'charlson',
                                  icd = 'icd10',
                                  assign0 = FALSE)
CCI <- score %>%
  select(PERSON_ID, wscore) %>% 
  dplyr::rename(cci = wscore) %>%
  mutate(cci_cat = ifelse(cci == 0, 0,
                          ifelse(cci == 1, 1, 
                                 ifelse(cci == 2, 2, "3+"))))
CCI$PERSON_ID <- as.numeric(CCI$PERSON_ID)

# 여기서, 첫검진년도 1년 전까지의 t40 자료가 없는 환자들은 0으로 부여 
df <- df %>% left_join(CCI)
df <- df %>% mutate(CCI = ifelse(is.na(cci), 0, cci),
                    CCI_CAT = ifelse(is.na(cci_cat), 0, cci_cat))

# CCI score 기록
score

rm(t40_CCI)
table(df$CCI_CAT)

save(list=ls(), file = "E:/seo_woman/age추가_fin/newbaseline_jk_gj_t206040.RData")

################################################################
## cox 분석을 위한 duration 추가
# 최초 날짜를 찾기위해 t20, t30, t40, t60 다 불러오기
rm(list=ls())
load("E:/seo_woman/age추가_fin/newbaseline_jk_gj_t206040.RData")

setwd('E:\\edata\\medicaldata\\30t')
listfile <- list.files()
t30 <- list()
for(i in 1:12){t30[[i]] <- data.table::fread(listfile[i],
                                             stringsAsFactors = F,
                                             encoding = "Latin-1") %>% 
  dplyr::select(KEY_SEQ, SEQ_NO,
                RECU_FR_DT, CLAUSE_CD
  ) %>% 
  dplyr::filter(KEY_SEQ %in% t20$KEY_SEQ)}
t30 <- dplyr::bind_rows(t30)
t30 <- t30 %>% left_join(t20 %>% select(PERSON_ID, KEY_SEQ))
t60 <- t60 %>% left_join(t20 %>% select(PERSON_ID, KEY_SEQ))

t20_first <- t20 %>% group_by(PERSON_ID) %>% summarise(fd = min(RECU_FR_DT))
t30_first <- t30 %>% group_by(PERSON_ID) %>% summarise(fd = min(RECU_FR_DT))
t40_first <- t40 %>% group_by(PERSON_ID) %>% summarise(fd = min(RECU_FR_DT))
t60_first <- t60 %>% group_by(PERSON_ID) %>% summarise(fd = min(RECU_FR_DT))

fd_t <- bind_rows(t20_first,t30_first,t40_first,t60_first) %>% 
  group_by(PERSON_ID) %>% summarise(firstdate = min(fd))

df <- df %>% left_join(fd_t)
df_f <- df %>% filter(is.na(firstdate))

## cox 분석을 위한 duration 구하기
## 사망환자: DTH_DATE - INDEX_DATE, 생존환자: 20131231 - INDEX_DATE
## 마지막으로 내원한 날짜 (final_date) 불러오기
t20 <- t20 %>% filter(PERSON_ID %in%df$PERSON_ID)
t20_last <- t20 %>% group_by(PERSON_ID) %>% summarise(FINAL_DATE = max(RECU_FR_DT))

df <- df %>% left_join(t20_last)
df$DTH_YM <- NULL

# 사망한 환자들 정보를 불러와서,
jk_DTH <- jk %>% filter(is.na(DTH_YM) == F)
df.dth <- df %>% filter(PERSON_ID %in% jk_DTH$PERSON_ID)
df.dth <- df.dth %>% left_join(jk_DTH %>% dplyr::select(PERSON_ID, DTH_YM))
# 마지막 내원날짜와 사망년월이 동일하면 마지막 내원 날짜(final_date)로, 아니면 사망년월(DTH_YM) 에 15일을 추가해 사망일자 생성
df.dth <- df.dth %>% mutate(DTH_DATE = ifelse(substr(FINAL_DATE, 1, 6) == DTH_YM,
                                              FINAL_DATE, as.integer(paste0(DTH_YM,15))))
df <- df %>% left_join(df.dth %>% select(PERSON_ID, DTH_DATE))

## duration 
# 1) INDEX_CANCER 2) DTH_DATE 3) 20131231
str(df)
df1 <- df %>% filter(CANCER == 1) %>% mutate(duration = as.integer(ymd(CANCER_DATE) - ymd(firstdate)))
df2 <- df %>% filter(CANCER == 0) %>% mutate(duration = ifelse(is.na(DTH_DATE) == F, as.integer(ymd(DTH_DATE) - ymd(firstdate)),
                                                               as.integer(ymd(20131231) - ymd(firstdate))))
df <- bind_rows(df1,df2)

save(list=ls(), file = "E:/seo_woman/age추가_fin/newbaseline_jk_gj_t20604030.RData")

###############################################################
##### 약제 정의 ########################
rm(list=ls())
load("E:/seo_woman/age추가_fin/newbaseline_jk_gj_t20604030.RData")

### PPI 정의
# firstyear 로 부터 3년 간 처방 받은 기록이 있으면 1, 아니면 0
ppi_drug <- c("2044", "1813", "6219", "4594", "2088", "2222", "5055")

t60_ppi <- t60 %>% 
  dplyr::mutate(ATC_CAT = substr(GNL_NM_CD, 1, 4)) %>% 
  dplyr::filter(ATC_CAT %in% ppi_drug) %>% 
  dplyr::left_join(t20 %>% 
                     dplyr::select(PERSON_ID,
                                   KEY_SEQ)) %>% 
  dplyr::left_join(df %>% dplyr::select(PERSON_ID, firstyear))

t60_ppi <- t60_ppi %>% filter(year(ymd(RECU_FR_DT)) <= firstyear + 2) # 3년 내에 기록이 있으며 ever 아니면 never
length(unique(t60_ppi$PERSON_ID)) 

df <- df %>% mutate(ppi = ifelse(PERSON_ID %in% t60_ppi$PERSON_ID, 1, 0))
table(df$ppi)

#### 진통제 정의
library(readxl)
medicine_p <- read_excel('E:\\edata\\seoprof\\newmetformin\\nsaid_code.xlsx')
medicine_p <- medicine_p %>% dplyr::select(Code,name)
medicine_p <- medicine_p %>% 
  dplyr::mutate(ATC_CAT = substr(Code, 1, 4)) %>% 
  unique() 

### 1단계 진통제 전체
t60_n <- t60 %>% 
  dplyr::mutate(ATC_CAT = substr(GNL_NM_CD, 1, 4)) %>% 
  dplyr::filter(ATC_CAT %in% medicine_p$ATC_CAT) %>% 
  dplyr::inner_join(medicine_p)
nsaid <- t60_n %>% 
  inner_join(t20 %>% dplyr::select(PERSON_ID, KEY_SEQ, RECU_FR_DT))  %>% 
  dplyr::left_join(df %>% dplyr::select(PERSON_ID, firstyear, CANCER_DATE)) %>% 
  arrange(PERSON_ID)

# 2년 단위를 뽑아서, 처방 빈도가 20번 이상이면 regular, 아니면 non-regular, 기록이 없으면 none
nsaid <- nsaid %>% filter(year(ymd(RECU_FR_DT)) - firstyear < 2)

nsaid <- nsaid %>% dplyr::select(-c(Code)) #%>% unique() %>% arrange(PERSON_ID, RECU_FR_DT)
nsaid_n <- nsaid %>% group_by(PERSON_ID, RECU_FR_DT, name) %>% summarise(kind_num = n(), day_sum = sum(MDCN_EXEC_FREQ))
nsaid_n <- nsaid_n %>% mutate(r_day = as.integer(day_sum / kind_num),
                              yearmon = substr(RECU_FR_DT,1,6))

nsaid_p_sum <- nsaid_n %>% group_by(PERSON_ID, yearmon) %>% summarise(t_1 = sum(day_sum), t_2 = sum(r_day))
nsaid_p_sum_1 <- nsaid_p_sum %>% mutate(month_t = ifelse(t_1>=15,1,0)) # 한달 단위 처방 일수
nsaid_p_sum_2 <- nsaid_p_sum %>% mutate(month_t = ifelse(t_2>=15,1,0)) # 한달 단위 처방 일수

# 진통약제를 15일 이상 받은 달이 6번 이상이면 환자면 regular?
nsaid_1 <- nsaid_p_sum_1 %>% group_by(PERSON_ID) %>% summarise(month = sum(month_t)) # 한달 단위 처방 일수
nsaid_1 <- nsaid_1 %>% mutate(n_rg = ifelse(month >= 6, 2, 1))
nsaid_2 <- nsaid_p_sum_2 %>% group_by(PERSON_ID) %>% summarise(month = sum(month_t)) # 한달 단위 처방 일수
nsaid_2 <- nsaid_2 %>% mutate(n_rg_2 = ifelse(month >= 3, 2, 1))

df <- df %>% left_join(nsaid_1 %>% select(PERSON_ID, n_rg))
df <- df %>% left_join(nsaid_2 %>% select(PERSON_ID, n_rg_2))

df <- df %>% mutate(nsaid = ifelse(is.na(n_rg), 0, n_rg),
                    nsaid_2 = ifelse(is.na(n_rg_2), 0, n_rg_2))


###########################
## 연속성도 한 번 봐보자 ##
###########################
nsaid_p_sum <- nsaid_n %>% group_by(PERSON_ID, yearmon) %>% summarise(t_1 = sum(day_sum), t_2 = sum(r_day))
nsaid_p_sum <- nsaid_p_sum %>% filter(t_2>=15) # 한달 단위 처방 일수

nsaid_p_list <- split(nsaid_p_sum, nsaid_p_sum$PERSON_ID)

# 한달에 15일 이상 처방을 n개월 연속으로 받은 환자 추출
case <- function(data,n) {
  idlist <- c()
  
  for(i in 1:length(data)){
    id <- data[[i]]
    list <- c()
    for(j in 1:(nrow(id) - n)) {
      list[j] <- as.numeric(data[[i]]$yearmon[j+n]) - as.numeric(data[[i]]$yearmon[j])
      if(sum(list == n, na.rm = T) >= 1){
        idlist[[i]] <- id
      }
    }
  }
  return(idlist)
}

aa <- case(nsaid_p_list,2)
aa <- bind_rows(aa)

length(unique(aa$PERSON_ID))

n_3 <- unique(aa$PERSON_ID) # 5315 명
nsaid_person <- unique(nsaid$PERSON_ID) # 5315 명

df <- df %>% mutate(nsaid3 = ifelse(PERSON_ID %in% n_3, 2, NA))
df <- df %>% mutate(nsaid_3 = ifelse(!PERSON_ID %in% nsaid_person, 0, 
                                     ifelse(is.na(nsaid3), 1, nsaid3)))

table(df$nsaid)
table(df$nsaid_2)
table(df$nsaid_3)


save(list=ls(), file = "E:/seo_woman/age추가_fin/final_dataframe_0207.RData")



#################################################################
################# 분석 진행 ##################################
rm(list=ls())
load("E:/seo_woman/age추가_fin/final_dataframe_0207.RData")
# load("E:/seo_woman/fin/newbaseline_fin_0202.RData")

# 콜레스테롤, 생존, 10단위 나이대 변수 생성
df <- df %>% mutate(tot_c = ifelse(TOT_CHOLE<240,0,1))
df <- df %>% mutate(survival = ifelse(is.na(DTH_DATE), 1, 0))
df <- df %>% mutate(age_b = ifelse(AGE_GROUP %in% c(6), 0,
                                   ifelse(AGE_GROUP %in% c(7,8), 1,
                                          ifelse(AGE_GROUP %in% c(9,10), 2,
                                                 ifelse(AGE_GROUP %in% c(11,12), 3, 4)))))
table(df$survival)
table(df$age_b)

# 연구대상자에서 의료급여, medicaid 해당 환자들 제외
df <- df %>% filter(income != 0 & insurance != 2)

## T1_new
df$CANCER <- factor(df$CANCER, levels = c(1,0))
table(df$CANCER)

xtabs(~df$AGE_GROUP+df$CANCER)
xtabs(~df$age_b+df$CANCER)
xtabs(~df$IPSN_TYPE_CD+df$CANCER)
xtabs(~df$income+df$CANCER)
xtabs(~df$bmi+df$CANCER, addNA =T)
xtabs(~df$tot_c+df$CANCER, addNA = T)
xtabs(~df$alcohol+df$CANCER, addNA = T)
#xtabs(~df$smoke+df$CANCER, addNA = T)
xtabs(~df$exer+df$CANCER, addNA = T)
xtabs(~df$hypertension+df$CANCER)
xtabs(~df$diabetes+df$CANCER)
xtabs(~df$CCI_CAT+df$CANCER)
# xtabs(~df$childbirth+df$CANCER)

xtabs(~df$nsaid+df$CANCER)
xtabs(~df$survival+df$CANCER)

# xtabs(~df$ace+df$CANCER)
# xtabs(~df$ibu+df$CANCER)
# xtabs(~df$dex+df$CANCER)


###### T2. nsaids 테이블 빈도 구하기 ######
### 1단계 진통제 전체
t60_n <- t60 %>% 
  dplyr::mutate(ATC_CAT = substr(GNL_NM_CD, 1, 4)) %>% 
  dplyr::filter(ATC_CAT %in% medicine_p$ATC_CAT) %>% 
  dplyr::inner_join(medicine_p)
nsaid <- t60_n %>% 
  inner_join(t20 %>% dplyr::select(PERSON_ID, KEY_SEQ, RECU_FR_DT))  %>% 
  dplyr::left_join(df %>% dplyr::select(PERSON_ID, firstyear, CANCER_DATE)) %>% 
  arrange(PERSON_ID)

# 2년 단위를 뽑아서, 
nd <- nsaid %>% filter(year(ymd(RECU_FR_DT)) - firstyear < 2)

## 1단계
acetaminophen <- nd %>% filter(grepl('acetaminophen', name))
length(unique(acetaminophen$PERSON_ID))
aspirin  <- nd %>% filter(grepl('aspirin', name)) 
length(unique(aspirin$PERSON_ID))
piroxicam  <- nd %>% filter(grepl('piroxicam', name)) 
length(unique(piroxicam$PERSON_ID))
diclofenac  <- nd %>% filter(grepl('diclofenac', name)) 
length(unique(diclofenac$PERSON_ID))
propacetamol  <- nd %>% filter(grepl('propacetamol', name)) 
length(unique(propacetamol$PERSON_ID))
etodolac  <- nd %>% filter(grepl('etodolac', name)) 
length(unique(etodolac$PERSON_ID))
celecoxib <- nd %>% filter(grepl('celecoxib', name)) 
length(unique(celecoxib$PERSON_ID))
ibuprofen <- nd %>% filter(grepl('ibuprofen', name)) 
length(unique(ibuprofen$PERSON_ID))
fenoprofen  <- nd %>% filter(grepl('fenoprofen', name)) 
length(unique(fenoprofen$PERSON_ID))
naproxen  <- nd %>% filter(grepl('naproxen', name)) 
length(unique(naproxen$PERSON_ID))
mefenamic  <- nd %>% filter(grepl('mefenamic', name)) 
length(unique(mefenamic$PERSON_ID))
nabumetone  <- nd %>% filter(grepl('nabumetone', name)) 
length(unique(nabumetone$PERSON_ID))
oxaprozin  <- nd %>% filter(grepl('oxaprozin', name)) 
length(unique(oxaprozin$PERSON_ID))
flurbiprofen  <- nd %>% filter(grepl('flurbiprofen', name)) 
length(unique(flurbiprofen$PERSON_ID))
ketorolac <- nd %>% filter(grepl('ketorolac', name)) 
length(unique(ketorolac$PERSON_ID))
clofenac <- nd %>% filter(grepl('clofenac', name)) 
length(unique(clofenac$PERSON_ID))
ketoprofen <- nd %>% filter(grepl('ketoprofen', name)) 
length(unique(ketoprofen$PERSON_ID))
rofecoxib <- nd %>% filter(grepl('rofecoxib', name)) 
length(unique(rofecoxib$PERSON_ID))
indomethacin <- nd %>% filter(grepl('indomethacin', name)) 
length(unique(indomethacin$PERSON_ID))
dexibuprofen <- nd %>% filter(grepl('dexibuprofen', name)) 
length(unique(dexibuprofen$PERSON_ID))


### T3. CCI table 채우기
str(score)
score <- score %>% filter(PERSON_ID %in% df$PERSON_ID)

sum(score$ami)
sum(score$chf)
sum(score$pvd)
sum(score$cevd)
sum(score$dementia)
sum(score$copd)
sum(score$rheumd)
sum(score$pud)
sum(score$mld)
sum(score$diab)
sum(score$diabwc)
sum(score$hp)
sum(score$rend)
sum(score$canc)
sum(score$msld)
sum(score$metacanc)
sum(score$aids)

## T4 chisq-test
# 분석을 위해 검진변수 NA 값 제거
df_d <- df[complete.cases(df[ , c("bmi", "alcohol", "smoke", "exer","tot_c")]), ]
table(df_d$IPSN_TYPE_CD)

# 일하는 지 여부(work) : 지역세대주, 직장가입자 에 해당하면 일함
df_d <- df_d %>% mutate(work = ifelse(IPSN_TYPE_CD %in% c(1,5), 1, 0))
df$CANCER <- factor(df$CANCER, levels = c(1,0))
table(df_d$work)

# 빈도 구하고,
table(df_d$CANCER)

xtabs(~age_b+CANCER, data = df_d)
xtabs(~AGE+CANCER, data = df_d)
xtabs(~IPSN_TYPE_CD+CANCER, data = df_d)
xtabs(~work+CANCER, data = df_d)
xtabs(~income+CANCER, data = df_d)

xtabs(~bmi+CANCER, data = df_d)
# xtabs(~tot_c+CANCER, data = df_d)
xtabs(~alcohol+CANCER, data = df_d)
# xtabs(~smoke+CANCER, data = df_d)
xtabs(~exer+CANCER, data = df_d)

# xtabs(~hypertension+CANCER, data = df_d)
# xtabs(~diabetes+CANCER, data = df_d)
xtabs(~CCI_CAT+CANCER, data = df_d)
# xtabs(~childbirth+CANCER, data = df_d)

xtabs(~nsaid+CANCER, data = df_d)
xtabs(~survival+CANCER, data = df_d)

# xtabs(~ace+CANCER, data = df_d)
# xtabs(~ibu+CANCER, data = df_d)
# xtabs(~dex+CANCER, data = df_d)

# chisq-test
chisq.test(xtabs(~age_b+CANCER, data = df_d))
chisq.test(xtabs(~IPSN_TYPE_CD+CANCER, data = df_d))
#chisq.test(xtabs(~work+CANCER, data = df_d))
chisq.test(xtabs(~income+CANCER, data = df_d))

chisq.test(xtabs(~bmi+CANCER, data = df_d))
# chisq.test(xtabs(~tot_c+CANCER, data = df_d))
chisq.test(xtabs(~alcohol+CANCER, data = df_d))
# chisq.test(xtabs(~smoke+CANCER, data = df_d))
chisq.test(xtabs(~exer+CANCER, data = df_d))

# chisq.test(xtabs(~hypertension+CANCER, data = df_d))
# chisq.test(xtabs(~diabetes+CANCER, data = df_d))
chisq.test(xtabs(~CCI_CAT+CANCER, data = df_d))
# chisq.test(xtabs(~childbirth+CANCER, data = df_d))

chisq.test(xtabs(~nsaid+CANCER, data = df_d))
chisq.test(xtabs(~survival+CANCER, data = df_d))


### T5_logistic
# 분석을 위해 검진변수 NA 값 제거

# 분석을 위한 factor 변환
df_d$CANCER <- factor(df_d$CANCER, levels = c(0,1))
df_d$AGE <- factor(df_d$AGE, levels = c(0,1,2))
# df_d$insurance <- factor(df_d$insurance, levels = c(0,1,2))
df_d$work <- factor(df_d$work, levels = c(0,1))
# df_d$IPSN_TYPE_CD <- factor(df_d$IPSN_TYPE_CD, levels = c(1,2,5,6,7,8))
# df_d$IPSN_TYPE_CD <- factor(df_d$IPSN_TYPE_CD, levels = c(2,1,6,5))
df_d$income <- factor(df_d$income, levels = c(1,2,3))
df_d$bmi <- factor(df_d$bmi, levels = c(0,1))
# df_d$tot_c <- factor(df_d$tot_c, levels = c(0,1))
df_d$alcohol <- factor(df_d$alcohol, levels = c(0,1))
#df_d$smoke <- factor(df_d$smoke, levels = c(0,1))
df_d$exer <- factor(df_d$exer, levels = c(0,1))
# df_d$hypertension <- factor(df_d$hypertension, levels = c(0,1))
# df_d$diabetes <- factor(df_d$diabetes, levels = c(0,1))
df_d$CCI_CAT <- factor(df_d$CCI_CAT, levels = c('0','1','2','3+'))
# df_d$childbirth <- factor(df_d$childbirth, levels = c(0,1))
df_d$nsaid <- factor(df_d$nsaid, levels = c(0,1,2))
df_d$nsaid_2 <- factor(df_d$nsaid_2, levels = c(0,1,2))
df_d$nsaid_3 <- factor(df_d$nsaid_3, levels = c(0,1,2))
df_d$survival <- factor(df_d$survival, levels = c(0,1))

# reference check
table(df_d$CANCER)
table(df_d$AGE)
table(df_d$work)
table(df_d$income)
mean(df_d$BMI)
table(df_d$alcohol)
table(df_d$exer)
table(df_d$CCI_CAT)
table(df_d$nsaid)
table(df_d$nsaid_2)
table(df_d$nsaid_3)
table(df_d$survival)

### logistic regression
fit <- glm(CANCER ~ AGE  + work +  income +
             bmi + alcohol + exer +
             CCI_CAT +
             nsaid+ survival, data=df_d, family = binomial(link="logit"))
summary(fit)
m1 <- tidy(fit, conf.int=T, exponentiate = T)
m1
summary(m1)
write.csv(m1, "E:/seo_woman/age추가_fin/logit_woman_fin_0208.csv")

### cox regression
table(df_d$CANCER)
table(df_d$work)

model1 <- coxph(Surv(duration, CANCER == 1) ~ work + AGE + income + 
                  BMI + alcohol + exer + 
                  CCI_CAT + nsaid + survival, data = df_d)
summary(model1)
write.csv(cbind(summary(model1)[[8]][,c(1,3,4)], summary(model1)[[7]][,5]),
          "E:/seo_woman/age추가_fin/cox_woman_fin.csv")

# nonworker, worker 군 나눠서
df_no <- df_d %>% filter(work == 0)
df_yes <- df_d %>% filter(work == 1)

model2 <- coxph(Surv(duration, CANCER == 1) ~ AGE + income + 
                  BMI + alcohol + exer + 
                  CCI_CAT + nsaid, data = df_no)
summary(model2)
write.csv(cbind(summary(model2)[[8]][,c(1,3,4)], summary(model2)[[7]][,5]),
          "E:/seo_woman/age추가_fin/cox_woman_no_fin.csv")

model2 <- coxph(Surv(duration, CANCER == 1) ~ AGE + income + 
                  BMI + alcohol + exer + 
                  CCI_CAT + nsaid, data = df_yes)
summary(model2)
write.csv(cbind(summary(model2)[[8]][,c(1,3,4)], summary(model2)[[7]][,5]),
          "E:/seo_woman/age추가_fin/cox_woman_yes_fin.csv")


### 생존 곡선 ###
df_no_c <- df_no %>% filter(CANCER == 1)
df_yes_c <- df_yes %>% filter(CANCER == 1)


ggsurvplot(survfit(Surv(duration/365, CANCER == 1) ~ nsaid, data=df_yes), data = df_yes, pval=FALSE, conf.int = TRUE, xlab = "Year",
           legend.labs = c("never", "non-regular", "regular"), legend.title = "Analgesics", linetype = "strata",
           risk.table = TRUE, risk.table.col = "strata", fun = "cumhaz", vpval.coord = c(0.1, 1)) 

ggsurvplot(survfit(Surv(duration/365, CANCER == 1) ~ nsaid, data=df_yes), data = df_yes, pval=FALSE, conf.int = TRUE, xlab = "Year",
           legend.labs = c("never", "non-regular", "regular"), legend.title = "Analgesics", linetype = "strata",
           risk.table = TRUE, risk.table.col = "strata", fun = "cumhaz",vpval.coord = c(0.1, 1))

ggrisktable(survfit(Surv(duration/365, CANCER == 1) ~ nsaid, data=df_no), data = df_no,
            survtable = "risk.table")


ee <- survfit(Surv(duration/365, CANCER == 1) ~ nsaid, data=df_no)
str(ee)
summary(ee)
nsaid_t <- data.frame(ee$n.risk, ee$cumhaz, ee$time)
colnames(nsaid_t) <- c("n.risk", "cumhaz", 'time')

nsaid_t
nw_ch0 <- nsaid_t %>% 
  filter(n.risk %in% c(1083, 1039, 944, 772, 14))
nw_ch1 <- nsaid_t %>% 
  filter(n.risk %in% c(3402,3345,	3184,	2955,	161))
nw_ch2 <- nsaid_t %>% 
  filter(n.risk %in% c(887,	878,833,765,97))


ee <- survfit(Surv(duration/365, CANCER == 1) ~ nsaid, data=df_yes)
str(ee)
summary(ee)
nsaid_t <- data.frame(ee$n.risk, ee$cumhaz, ee$time)
colnames(nsaid_t) <- c("n.risk", "cumhaz", 'time')

nsaid_t
nw_ch0 <- nsaid_t %>% 
  filter(n.risk %in% c(432,405,338,261,5))
nw_ch1 <- nsaid_t %>% 
  filter(n.risk %in% c(1054,1031,992,912,37))
nw_ch2 <- nsaid_t %>% 
  filter(n.risk %in% c(229,222,211,197,23))

cancer <- df_d %>% filter(CANCER == 1)
xtabs(~nsaid + AGE, data = cancer)
ctest <- chisq.test(xtabs(~nsaid + AGE, data = cancer))
ctest
ctest$stdres
ctest$expecte

control <- df_d %>% filter(CANCER == 0)
xtabs(~nsaid + AGE, data = control)
ctest <- chisq.test(xtabs(~nsaid + AGE, data = control))
ctest
ctest$stdres
ctest$expecte

xtabs(~nsaid + AGE, data = df_d)
ctest <- chisq.test(xtabs(~nsaid + AGE, data = df_d))
ctest
ctest$stdres
ctest$expecte

# hazard ratio plot
library(survival)
library(scales)
# if(!require(devtools)) install.packages("devtools")
# devtools::install_github("kassambara/survminer", build_vignettes = TRUE)# install.packages("survminer")
library("survminer")

fit2 <- coxph(Surv(duration, CANCER==1) ~ AGE + income + 
                BMI + alcohol + exer + 
                CCI_CAT + nsaid + survival, data = df_no)
ggforest(fit2) 

fit2 <- coxph(Surv(duration, CANCER==1) ~ AGE + income + 
                BMI + alcohol + exer + 
                CCI_CAT + nsaid + survival, data = df_yes)
ggforest(fit2) 


### work, non-work 군으로 나눠서 빈도
df_d$CANCER <- factor(df_d$CANCER, levels = c(1,0))
xtabs(~CANCER+work, data = df_d)
xtabs(~AGE+work, data = df_d)
xtabs(~income+work, data = df_d)

xtabs(~bmi+work, data = df_d)
xtabs(~alcohol+work, data = df_d)
xtabs(~exer+work, data = df_d)

xtabs(~CCI_CAT+work, data = df_d)

xtabs(~nsaid+work, data = df_d)
xtabs(~survival+work, data = df_d)


#####################
### age 별로 분석 ###
#####################
table(df_d$CANCER)

df_d$CANCER <- factor(df_d$CANCER, levels = c(0,1))
df_1 <- df_d %>% filter(AGE %in% c(0))
df_2 <- df_d %>% filter(AGE %in% c(1))
df_3 <- df_d %>% filter(AGE %in% c(2))

table(df_1$CANCER)
table(df_1$work)
table(df_1$income)
mean(df_1$BMI)
table(df_1$alcohol)
table(df_1$exer)
table(df_1$CCI_CAT)
table(df_1$nsaid)

# logit
fit_1 <- glm(CANCER ~ work +  income +
               bmi + alcohol + exer +
               CCI_CAT +
               nsaid, data=df_1, family = binomial(link="logit"))
summary(fit_1)
m1 <- tidy(fit_1, conf.int=T, exponentiate = T)
m1
summary(m1)
write.csv(m1, "E:/seo_woman/age추가_fin/logit_woman_1_fin.csv")

# cox
model1 <- coxph(Surv(duration, CANCER == 1) ~ work + income + 
                  BMI + alcohol + exer + 
                  CCI_CAT + nsaid, data = df_1)
summary(model1)
write.csv(cbind(summary(model1)[[8]][,c(1,3,4)], summary(model1)[[7]][,5]),
          "E:/seo_woman/age추가_fin/cox_woman_age1_fin.csv")

# 30 nonworker, worker 
df_1_no <- df_1 %>% filter(work == 0)
df_1_yes <- df_1 %>% filter(work == 1)

model2 <- coxph(Surv(duration, CANCER == 1) ~ income + 
                  BMI + alcohol + exer + 
                  CCI_CAT + nsaid, data = df_1_no)
summary(model2)
write.csv(cbind(summary(model2)[[8]][,c(1,3,4)], summary(model2)[[7]][,5]),
          "E:/seo_woman/age추가_fin/cox_woman_age1_no_fin.csv")

model2 <- coxph(Surv(duration, CANCER == 1) ~ income + 
                  BMI + alcohol + exer + 
                  CCI_CAT + nsaid, data = df_1_yes)
summary(model2)
write.csv(cbind(summary(model2)[[8]][,c(1,3,4)], summary(model2)[[7]][,5]),
          "E:/seo_woman/age추가_fin/cox_woman_age1_yes_fin.csv")

df_1$CANCER <- factor(df_1$CANCER, levels = c(1,0))
table(df_1$CANCER)

xtabs(~work+CANCER, data = df_1)
xtabs(~income+CANCER, data = df_1)

xtabs(~bmi+CANCER, data = df_1)
xtabs(~alcohol+CANCER, data = df_1)
xtabs(~exer+CANCER, data = df_1)

xtabs(~CCI_CAT+CANCER, data = df_1)

xtabs(~nsaid+CANCER, data = df_1)
xtabs(~survival+CANCER, data = df_1)

# work 로 나눠서 빈도 구하기
table(df_1$work)
xtabs(~CANCER+work, data = df_1)
xtabs(~income+work, data = df_1)

xtabs(~bmi+work, data = df_1)
xtabs(~alcohol+work, data = df_1)
xtabs(~exer+work, data = df_1)

xtabs(~CCI_CAT+work, data = df_1)

xtabs(~nsaid+work, data = df_1)
xtabs(~survival+work, data = df_1)

## 40대 환자
table(df_2$CANCER)
table(df_2$work)
table(df_2$income)
mean(df_2$BMI)
table(df_2$alcohol)
table(df_2$exer)
table(df_2$CCI_CAT)
table(df_2$nsaid)

#logit
fit_2 <- glm(CANCER ~ work +  income +
               bmi + alcohol + exer +
               CCI_CAT +
               nsaid + survival, data=df_2, family = binomial(link="logit"))
summary(fit_2)
m2 <- tidy(fit_2, conf.int=T, exponentiate = T)
m2
summary(m2)
write.csv(m2, "E:/seo_woman/age추가_fin/logit_woman_2_fin.csv")

# cox
model1 <- coxph(Surv(duration, CANCER == 1) ~ work + income + 
                  BMI + alcohol + exer + 
                  CCI_CAT + nsaid, data = df_2)
summary(model1)
write.csv(cbind(summary(model1)[[8]][,c(1,3,4)], summary(model1)[[7]][,5]),
          "E:/seo_woman/age추가_fin/cox_woman_age2_fin.csv")

# 40 nonworker, worker 
df_2_no <- df_2 %>% filter(work == 0)
df_2_yes <- df_2 %>% filter(work == 1)

model2 <- coxph(Surv(duration, CANCER == 1) ~ income + 
                  BMI + alcohol + exer + 
                  CCI_CAT + nsaid, data = df_2_no)
summary(model2)
write.csv(cbind(summary(model2)[[8]][,c(1,3,4)], summary(model2)[[7]][,5]),
          "E:/seo_woman/age추가_fin/cox_woman_age2_no_fin.csv")

model2 <- coxph(Surv(duration, CANCER == 1) ~ income + 
                  BMI + alcohol + exer + 
                  CCI_CAT + nsaid, data = df_2_yes)
summary(model2)
write.csv(cbind(summary(model2)[[8]][,c(1,3,4)], summary(model2)[[7]][,5]),
          "E:/seo_woman/age추가_fin/cox_woman_age2_yes_fin.csv")

df_2$CANCER <- factor(df_2$CANCER, levels = c(1,0))
table(df_2$CANCER)

xtabs(~work+CANCER, data = df_2)
xtabs(~income+CANCER, data = df_2)

xtabs(~bmi+CANCER, data = df_2)
xtabs(~alcohol+CANCER, data = df_2)
xtabs(~exer+CANCER, data = df_2)

xtabs(~CCI_CAT+CANCER, data = df_2)

xtabs(~nsaid+CANCER, data = df_2)
xtabs(~survival+CANCER, data = df_2)


## 50대 환자
table(df_3$CANCER)
table(df_3$work)
table(df_3$income)
mean(df_3$BMI)
table(df_3$alcohol)
table(df_3$exer)
table(df_3$CCI_CAT)
table(df_3$nsaid)

fit_2 <- glm(CANCER ~ work +  income +
               bmi + alcohol + exer +
               CCI_CAT +
               nsaid + survival, data=df_3, family = binomial(link="logit"))
summary(fit_2)
m2 <- tidy(fit_2, conf.int=T, exponentiate = T)
m2
summary(m2)
write.csv(m2, "E:/seo_woman/age추가_fin/logit_woman_3_fin.csv")

# cox
model1 <- coxph(Surv(duration, CANCER == 1) ~ work + income + 
                  BMI + alcohol + exer + 
                  CCI_CAT + nsaid, data = df_3)
summary(model1)
write.csv(cbind(summary(model1)[[8]][,c(1,3,4)], summary(model1)[[7]][,5]),
          "E:/seo_woman/age추가_fin/cox_woman_age3_fin.csv")

# 50 nonworker, worker 
df_3_no <- df_3 %>% filter(work == 0)
df_3_yes <- df_3 %>% filter(work == 1)

model2 <- coxph(Surv(duration, CANCER == 1) ~ income + 
                  BMI + alcohol + exer + 
                  CCI_CAT + nsaid, data = df_3_no)
summary(model2)
write.csv(cbind(summary(model2)[[8]][,c(1,3,4)], summary(model2)[[7]][,5]),
          "E:/seo_woman/age추가_fin/cox_woman_age3_no_fin.csv")

model2 <- coxph(Surv(duration, CANCER == 1) ~ income + 
                  BMI + alcohol + exer + 
                  CCI_CAT + nsaid, data = df_3_yes)
summary(model2)
write.csv(cbind(summary(model2)[[8]][,c(1,3,4)], summary(model2)[[7]][,5]),
          "E:/seo_woman/age추가_fin/cox_woman_age3_yes_fin.csv")

df_3$CANCER <- factor(df_3$CANCER, levels = c(1,0))
table(df_3$CANCER)

xtabs(~work+CANCER, data = df_3)
xtabs(~income+CANCER, data = df_3)

xtabs(~bmi+CANCER, data = df_3)
xtabs(~alcohol+CANCER, data = df_3)
xtabs(~exer+CANCER, data = df_3)

xtabs(~CCI_CAT+CANCER, data = df_3)

xtabs(~nsaid+CANCER, data = df_3)
xtabs(~survival+CANCER, data = df_3)




myd <- data.frame(   var1  = c(1,1,2,2,3,3,4,4),  
                     samp = c("A","B","A","B","A","B","A","B"),  
                     Value1 = c(3.5,2,5,8,3,2,7,2), Value2 = c(1.5,0,5,5,3,0,4,5) )
# rshaping data to long form for ggplot2
library(reshape2)
meltd<- melt(myd, id.vars=1:2)

#plot
library(ggplot2)
ggplot(meltd, aes(x=var1, y=value, fill=variable)) +
  geom_bar(stat="identity") + facet_grid(~samp) + theme_bw()



gp <- read.csv("E:/seo_woman/graph_stackedbar.csv")
gp$cancer <- factor(gp$cancer, levels = c('Case', 'Control', 'All'))
gp$nsaid <- factor(gp$nsaid, levels = c('Regular', 'Non-regular', 'Never'))

ggplot(gp, aes(x=cancer, y=ratio, fill=nsaid), size=4) +
  geom_bar(stat="identity", position = position_stack(reverse = T)) + facet_grid(~age) + theme_bw()  + 
  geom_text(aes(y = ratio, label = paste0(ratio,'%'), size = 3, family="serif"), 
            position = position_stack(vjust = 0.5, reverse = T), size = 4, color = "black") + ylab("") + xlab("") + 
  theme(plot.title = element_text(size = 18, family = "serif", face = "bold"),
        text=element_text(family="serif"),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10)) +
  guides(fill=guide_legend(title="Analgesics", reverse =T)) +
  theme(strip.text = element_text(size=12))+
  theme(axis.line = element_line(size=0.6, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank())  +
  scale_fill_hue(c=40, l=88) #+ theme(legend.title=element_blank())
 











