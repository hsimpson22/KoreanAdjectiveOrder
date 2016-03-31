#Korean Adjective order final paper

### I want to test whether the order of prenominal adjectives in Korean is determined, among other things, by their lengths. This is a script that
# - retrieves all instances of two adjectives preceding a noun
# - computes the lengths of the adjectives in both positions as well as the pairwise differences
# - plots the data and computes a significance test
#First, in bash shell convert KAIST files to UTF-8
#for i in *; do iconv -f EUC-KR -t UTF-8 $i > $i.new; done
#mkdir non_unicode/
#mv *.txt non_unicode/
#for i in *.new; do mv $i `basename $i .new`; done
# clear memory
rm(list=ls(all=TRUE))

#Loading script for nice regex searching with greedy return of string matches
source("http://www.linguistics.ucsb.edu/faculty/stgries/exact_matches.r")

# determine the corpus files
kaist.files <- dir(".", full.names=TRUE) #
kaist.mod.mod.n <-c()
kaist.file<-c()
for (i in 1:length(kaist.files)){
  cat(i, "\n")
  curr.kaist.file<-scan(kaist.files[i], what=character(0), sep="\n", quiet=TRUE)
  # remove metadata lines
  curr.kaist.file<-gsub("<.*", "", curr.kaist.file, perl=TRUE)
  curr.kaist.matches <-exact.matches("[^\\s]+\\s+([^\\s]+n[^\\s]+하/...(.)?+[^는].etm|[^\\s]+xsn(..)?|[^\\s]+paa([^\\s]+(etm|ecc))?)\\s+[^\\s]+\\s+([^\\s]+n[^\\s]+하/...(.)?+[^는].etm|[^\\s]+xsn(..)?|[^\\s]+paa([^\\s]+etm)?)\\s+[^\\s]+\\s+[^\\s]+n[^\\s]*", curr.kaist.file, pcre=TRUE, case.sens=FALSE)[[1]]
  kaist.file<-c(kaist.file, rep(kaist.files[i], length(curr.kaist.matches)))
  kaist.mod.mod.n <- c(kaist.mod.mod.n, curr.kaist.matches)
}  
length(exact.matches("\t있/paa", kaist.mod.mod.n, pcre=TRUE, case.sens=FALSE)[[1]])#45
length(exact.matches("\t아니/paa", kaist.mod.mod.n, pcre=TRUE, case.sens=FALSE)[[1]])#113
head(kaist.file); length(kaist.file) 
kaist.mod.mod.n2<-paste(kaist.file, kaist.mod.mod.n, sep="\t")
kaist.mod.mod.n2<-kaist.mod.mod.n2[(grepl("\t있/", kaist.mod.mod.n2))==FALSE]#so something lexicalized like 맛있 is still maintained
kaist.mod.mod.n2<-kaist.mod.mod.n2[(grepl("\t없/", kaist.mod.mod.n2))==FALSE]
kaist.mod.mod.n2<-kaist.mod.mod.n2[(grepl("\t아니/", kaist.mod.mod.n2))==FALSE]
#remove 같/paa+은/etm"?? can you say 같고 길은 사람
head(kaist.mod.mod.n2); length(kaist.mod.mod.n2) #3911
write.table(kaist.mod.mod.n2, file=("../KAIST_raw_pos_modmodn.csv"),quote=F, sep="\t")
kaist.analyzed <-read.table(file=("../KAIST_raw_pos_modmodn.csv"), quote="", sep="\t", header=TRUE, comment.char="")
str(kaist.analyzed)

#Add one more set of mod mod n matches with mma included (START)
 kaist.files <- dir(".", full.names=TRUE) #
 kaist.mma.matches <-c()
 kaist.file<-c()
for (i in 1:length(kaist.files)){
  cat(i, "\n")
  curr.kaist.file<-scan(kaist.files[i], what=character(0), sep="\n", quiet=TRUE)
  # remove metadata lines
  curr.kaist.file<-gsub("<.*", "", curr.kaist.file, perl=TRUE) 
  curr.kaist.file<-curr.kaist.file[which(grepl("mma", curr.kaist.file)==TRUE)]
  curr.mma.matches<-exact.matches("[^\\s]+\\s+([^\\s]+n[^\\s]+하/...(.)?+[^는].etm|[^\\s]+xsn(..)?|[^\\s]+mma|[^\\s]+paa([^\\s]+(etm|ecc))?)\\s+[^\\s]+\\s+([^\\s]+n[^\\s]+하/...(.)?+[^는].etm|[^\\s]+xsn(..)?|[^\\s]+mma|[^\\s]+paa([^\\s]+etm)?)\\s+[^\\s]+\\s+[^\\s]+n[^\\s]*", curr.kaist.file, pcre=TRUE, case.sens=FALSE)[[1]]
  kaist.file<-c(kaist.file, rep(kaist.files[i], length(curr.mma.matches)))
  kaist.mma.matches <- c(kaist.mma.matches, curr.mma.matches)  
}
length(kaist.mma.matches);head(kaist.mma.matches)
length(kaist.file)
kaist.mma.matches<-paste(kaist.file, kaist.mma.matches, sep="\t")
kaist.mma.matches<-(kaist.mma.matches[which(grepl("mma", kaist.mma.matches)==TRUE)])
kaist.mma.matches <-kaist.mma.matches[(grepl("\\t있/", kaist.mma.matches, perl=TRUE))==FALSE]
kaist.mma.matches <-kaist.mma.matches[(grepl("\\t없/", kaist.mma.matches, perl=TRUE))==FALSE]
kaist.mma.matches <-kaist.mma.matches[(grepl("\\t아니/", kaist.mma.matches, perl=TRUE))==FALSE]
#Add removal of "같" 
kaist.mma.matches <-kaist.mma.matches[(grepl("\\t같/", kaist.mma.matches, perl=TRUE))==FALSE] 
length(kaist.mma.matches); head(kaist.mma.matches)

write.table(kaist.mma.matches, file=("../KAIST_raw_pos_mma.csv"),quote=F, sep="\t")
kaist.mma<-read.table(file=("../KAIST_raw_pos_mma.csv"), quote="", sep="\t", header=TRUE, comment.char="")
#see if it is possible for mma to precede another type of modifier
table(grepl("mma", kaist.mma$MOD2_ANALYSIS[which(grepl("mma", kaist.mma$MOD1_ANALYSIS)==TRUE)])) #FALSE=647, TRUE = 667
MOD1.mma<-kaist.mma[which(grepl("mma", kaist.mma$MOD1_ANALYSIS)==TRUE),]
MOD1.mma<-MOD1.mma[which(grepl("mma", MOD1.mma$MOD2_ANALYSIS, perl=TRUE)==FALSE),]
MOD2.mma<-kaist.mma[which(grepl("mma", kaist.mma$MOD2_ANALYSIS)==TRUE),]
MOD2.mma<-MOD2.mma[which(grepl("mma", MOD2.mma$MOD1_ANALYSIS, perl=TRUE)==FALSE),]
length(union(unique(MOD1.mma$MOD1), unique(MOD2.mma$MOD2))) #45
mma.list<-(union(unique(MOD1.mma$MOD1), unique(MOD2.mma$MOD2)))
MMA.PHRASE<-(paste(kaist.mma$MOD1, kaist.mma$MOD2, kaist.mma$NOUN, sep=" "));head(MMA.PHRASE)
#Get rid of rid of weird punctuation ["()[]\}&,]
mma.list <-gsub("[\\\\\\,&}\\.\\(\\)\\\"\\[\\]]", "", mma.list, perl=TRUE)
mma.list <-unique(mma.list)
mma.ex<-c()
for (i in 1:(length(mma.list))) {
	mma.ex.curr<-MMA.PHRASE[which(grepl(mma.list[i], MMA.PHRASE))==TRUE]	
	mma.ex<-c(mma.ex, mma.ex.curr)
	}
write.table(mma.list, file=("../KAIST_mmalist.csv"),quote=F, sep="\t")
length(intersect(unique(MOD1.mma$MOD1), unique(MOD2.mma$MOD2)))
length(union(unique(MOD1.mma$MOD2), unique(MOD2.mma$MOD1)))
length(intersect(unique(MOD1.mma$MOD2), unique(MOD2.mma$MOD1)))
#Add one more set of mod mod n matches with mma included (END)

#Add one more set of mod mod n matches with noun + 의 constructions included (START)
kaist.files <- dir(".", full.names=TRUE) #
kaist.genitive.matches <-c()
kaist.file<-c()
for (i in 1:length(kaist.files)){
  cat(i, "\n")
  curr.kaist.file<-scan(kaist.files[i], what=character(0), sep="\n", quiet=TRUE)
  # remove metadata lines
  curr.kaist.file<-gsub("<.*", "", curr.kaist.file, perl=TRUE)
#adding plain noun version:  curr.genitive.matches <-exact.matches("[^\\s]+\\s+([^\\s]+n[^\\s]+하/...(.)?+[^는].etm|[^\\s]+xsn(..)?|[^\\s]+mma|[^\\s]+ncn(.의.jcm)?|[^\\s]+paa([^\\s]+(etm|ecc))?)\\s+[^\\s]+\\s+([^\\s]+n[^\\s]+하/...(.)?+[^는].etm|[^\\s]+xsn(..)?|[^\\s]+mma|[^\\s]+ncn(.의.jcm)?|[^\\s]+paa([^\\s]+etm)?)\\s+[^\\s]+\\s+[^\\s]+n[^\\s]*", curr.kaist.file, pcre=TRUE, case.sens=FALSE)[[1]]
	curr.genitive.matches <-exact.matches("[^\\s]+\\s+([^\\s]+n[^\\s]+하/...(.)?+[^는].etm|[^\\s]+xsn(..)?|[^\\s]+mma|[^\\s]+ncn.의.jcm|[^\\s]+paa([^\\s]+(etm|ecc))?)\\s+[^\\s]+\\s+([^\\s]+n[^\\s]+하/...(.)?+[^는].etm|[^\\s]+xsn(..)?|[^\\s]+mma|[^\\s]+ncn.의.jcm|[^\\s]+paa([^\\s]+etm)?)\\s+[^\\s]+\\s+[^\\s]+n[^\\s]*", curr.kaist.file, pcre=TRUE, case.sens=FALSE)[[1]]
  kaist.file<-c(kaist.file, rep(kaist.files[i], length(curr.genitive.matches)))
  kaist.genitive.matches <- c(kaist.genitive.matches, curr.genitive.matches)
}  
head(kaist.file); length(kaist.file) 
head(kaist.genitive.matches); length(kaist.genitive.matches) #56732
kaist.genitive.matches.orig<-kaist.genitive.matches
kaist.file.orig<-kaist.file
kaist.genitive.matches.orig<-paste(kaist.file.orig, kaist.genitive.matches.orig, sep="\t")
kaist.genitive.matches<-kaist.genitive.matches.orig
kaist.genitive.matches <-(kaist.genitive.matches[which(grepl("jcm\t", kaist.genitive.matches)==TRUE)]) #this should get non final noun + 의 (so must be one of the two modifiers)
head(kaist.genitive.matches); length(kaist.genitive.matches) #49495
kaist.genitive.matches<-kaist.genitive.matches[(grepl("\t있/", kaist.genitive.matches))==FALSE]
kaist.genitive.matches<-kaist.genitive.matches[(grepl("\t없/", kaist.genitive.matches))==FALSE]
kaist.genitive.matches<-kaist.genitive.matches[(grepl("\t아니/", kaist.genitive.matches))==FALSE]
kaist.genitive.matches <-kaist.genitive.matches[(grepl("\t같/", kaist.genitive.matches, perl=TRUE))==FALSE]
head(kaist.genitive.matches); length(kaist.genitive.matches) #45037
#kaist.noungen.matches<-kaist.genitive.matches
#kaist.genitive.matches<-(kaist.genitive.matches[which(grepl("의.jcm\t", kaist.genitive.matches)==TRUE)]) 
#head(kaist.genitive.matches); length(kaist.genitive.matches)  #186100
#get rid of plain nouns
#kaist.genitive.matches <-(kaist.genitive.matches[which(grepl("ncn\t", kaist.genitive.matches)==FALSE)]) 
#head(kaist.genitive.matches); length(kaist.genitive.matches)
write.table(kaist.genitive.matches, file=("../KAIST_raw_pos_genitive.csv"),quote=F, sep="\t")
kaist.gen <-read.table(file=("../KAIST_raw_pos_genitive.csv"), quote="", sep="\t", header=TRUE, comment.char="") 
str(kaist.gen) 
#see how often gen follows or precedes another type of modifier
table(grepl("의/jcm", kaist.gen$MOD2_ANALYSIS[which(grepl("의/jcm", kaist.gen$MOD1_ANALYSIS)==TRUE)])) #of triples with genitive constructions in MOD1 position, how many also have a genitive construction in MOD2 position: FALSE= 14571, TRUE = 11817
MOD1.gen<-kaist.gen[which(grepl("의/jcm", kaist.gen$MOD1_ANALYSIS)==TRUE),]
MOD1.gen<-MOD1.gen[which(grepl("의/jcm", MOD1.gen$MOD2_ANALYSIS, perl=TRUE)==FALSE),] #the genitive constructions in MOD1 position where MOD2 is NOT a genitive construction, 
str(MOD1.gen)#14571
MOD2.gen<-kaist.gen[which(grepl("의/jcm", kaist.gen$MOD2_ANALYSIS)==TRUE),]
MOD2.gen<-MOD2.gen[which(grepl("의/jcm", MOD2.gen$MOD1_ANALYSIS, perl=TRUE)==FALSE),]#the genitive constructions in MOD2 position where MOD1 is NOT a genitive construction
str(MOD2.gen) #18649
#find out how many modifier types there are that occur in either position
length(union(unique(MOD1.gen$MOD1), unique(MOD2.gen$MOD2))) #10934
#find out how many modifier types there are that occur in both positions
length(intersect(unique(MOD1.gen$MOD1), unique(MOD2.gen$MOD2)))
#find out how many modifier types there are that 
length(union(unique(MOD1.gen$MOD2), unique(MOD2.gen$MOD1)))
length(intersect(unique(MOD1.gen$MOD2), unique(MOD2.gen$MOD1)))
gen.versatile<-(intersect(unique(MOD1.gen$MOD1), unique(MOD2.gen$MOD2)))
write.table(gen.versatile, file=("../KAIST_GenVersatile_Review.csv"),quote=F, sep="\t")
#Add one more set of mod mod n matches with noun + 의 constructions included (START)

#CREATE CODING FILE FOR WONA (START)
#PHRASE<-(paste(MOD1, MOD2, NOUN, sep=" "));head(PHRASE)
#kaist.wona.coding<-data.frame(PHRASE, MOD1, MOD2,NOUN)
#write.table(kaist.wona.coding, file=("../KAIST_Wona_Coding.csv"),quote=F, sep="\t")
#CREATE FILE FOR WONA (END)

#CREATE CODING FILE FOR ME (START)
#PHRASE_ANALYSIS<-(paste(MOD1_ANALYSIS, MOD2_ANALYSIS, NOUN_ANALYSIS, sep=" "));head(PHRASE_ANALYSIS)
#kaist.me.coding<-data.frame(PHRASE[], MOD1, MOD2,NOUN)
#write.table(kaist.me.coding, file=("../KAIST_Me_Coding.csv"),quote=F, sep="\t")
#CREATE CODING FILE FOR ME (END)

#READ BACK IN CODING FILE (START)
#kaist.coded <-read.table(file=("../KAIST_Combined_Coding.csv"), quote="", sep="\t", header=TRUE, comment.char="")
kaist.coded <-read.table(file=("../KAIST_Wona_Coding_final.csv"), quote="", sep="\t", header=TRUE, comment.char="")
str(kaist.coded)
#attach both bc the coding columns are not in kaist.analyzed so won't overwrite with empty
attach(kaist.coded)
attach(kaist.analyzed)
#add analysis columns back in from KAIST raw mod mod n
kaist.mod.mod.n<-data.frame(GOOD, MOD1, MOD1_ANALYSIS, MOD1_DEF,MOD1_SUBJOBJ, MOD1_AFFLOAD, MOD2, MOD2_ANALYSIS, MOD2_DEF, MOD2_SUBJOBJ, MOD2_AFFLOAD, NOUN, NOUN_ANALYSIS,NOUN_DEF)
str(kaist.mod.mod.n)
detach(kaist.coded)
detach(kaist.analyzed)
good.matchn<-grepl("n", kaist.mod.mod.n$GOOD, perl=TRUE, ignore.case=TRUE)
kaist.mod.mod.n<-kaist.mod.mod.n[which(good.matchn=="FALSE"),]
library(gdata)
kaist.mod.mod.n<-drop.levels(kaist.mod.mod.n)
good.matchcheck<-grepl("check", kaist.mod.mod.n$GOOD, perl=TRUE, ignore.case=TRUE)
kaist.mod.mod.n<-kaist.mod.mod.n[which(good.matchcheck=="FALSE"),]
kaist.mod.mod.n<-drop.levels(kaist.mod.mod.n)
subjobj1.match<-grepl("\\d", kaist.mod.mod.n$MOD1_SUBJOBJ, perl=TRUE, ignore.case=TRUE) 
kaist.mod.mod.n<-kaist.mod.mod.n[which(subjobj1.match=="TRUE"),]
kaist.mod.mod.n<-drop.levels(kaist.mod.mod.n)
subjobj2.match<-grepl("\\d", kaist.mod.mod.n$MOD2_SUBJOBJ, perl=TRUE, ignore.case=TRUE)
kaist.mod.mod.n<-kaist.mod.mod.n[which(subjobj2.match=="TRUE"),]
kaist.mod.mod.n<-drop.levels(kaist.mod.mod.n)
affload1.match<-grepl("\\d", kaist.mod.mod.n$MOD1_AFFLOAD, perl=TRUE, ignore.case=TRUE)
kaist.mod.mod.n<-kaist.mod.mod.n[which(affload1.match=="TRUE"),]
kaist.mod.mod.n<-drop.levels(kaist.mod.mod.n)
affload2.match<-grepl("\\d", kaist.mod.mod.n$MOD2_AFFLOAD, perl=TRUE, ignore.case=TRUE)
kaist.mod.mod.n<-kaist.mod.mod.n[which(affload2.match=="TRUE"),]
kaist.mod.mod.n<-drop.levels(kaist.mod.mod.n)
str(kaist.mod.mod.n)

attach(kaist.mod.mod.n)
sum(grepl("\\(|\\)", MOD1, perl=TRUE))
sum(grepl("\\(|\\)", MOD2, perl=TRUE))
sum(grepl("\\(|\\)", NOUN, perl=TRUE))
#kaist.mod.mod.n<-kaist.mod.mod.n[grepl(kaist.mod.mod.n$GOOD)]
#clean up /sl (slang) tags 
MOD1_ANALYSIS <-gsub("(^|\t).?/sl\\+", "", MOD1_ANALYSIS, perl=TRUE)
MOD2_ANALYSIS <-gsub("(^|\t).?/sl\\+", "", MOD2_ANALYSIS, perl=TRUE)
NOUN_ANALYSIS <-gsub("(^|\t).?/sl\\+", "", NOUN_ANALYSIS, perl=TRUE)
#get rid of anything within parentheses in analysis columns (so roots can be compared)
MOD1_ANALYSIS <-gsub("\\(.*\\)", "", MOD1_ANALYSIS, perl=TRUE);length(MOD1_ANALYSIS)
MOD2_ANALYSIS <-gsub("\\(.*\\)", "", MOD2_ANALYSIS, perl=TRUE);length(MOD2_ANALYSIS)
NOUN_ANALYSIS <-gsub("\\(.*\\)", "", NOUN_ANALYSIS, perl=TRUE);length(NOUN_ANALYSIS)
#or on one side of parentheses
MOD1_ANALYSIS <-gsub("\\([^\\)]+$", "", MOD1_ANALYSIS, perl=TRUE);length(MOD1_ANALYSIS)
MOD2_ANALYSIS <-gsub("\\([^\\)]+$", "", MOD2_ANALYSIS, perl=TRUE);length(MOD2_ANALYSIS)
NOUN_ANALYSIS <-gsub("\\([^\\)]+$", "", NOUN_ANALYSIS, perl=TRUE);length(NOUN_ANALYSIS)
#Get rid of rid of weird punctuation ["()[]\}&,]
MOD1_ANALYSIS <-gsub("[\\\\\\,&}\\.\\(\\)\\\"\\[\\]]", "", MOD1_ANALYSIS, perl=TRUE)
MOD2_ANALYSIS <-gsub("[\\\\\\,&}\\.\\(\\)\\\"\\[\\]]", "", MOD2_ANALYSIS, perl=TRUE)
NOUN_ANALYSIS <-gsub("[\\\\\\,&}\\.\\(\\)\\\"\\[\\]]", "", NOUN_ANALYSIS, perl=TRUE)
#clean up /sl (slang) tags 
MOD1 <-gsub("(^|\t).?/sl\\+", "", MOD1, perl=TRUE)
MOD2 <-gsub("(^|\t).?/sl\\+", "", MOD2, perl=TRUE)
NOUN <-gsub("(^|\t).?/sl\\+", "", NOUN, perl=TRUE)
#get rid of anything within parentheses in analysis columns (so roots can be compared)
MOD1 <-gsub("\\(.*\\)", "", MOD1, perl=TRUE);length(MOD1)
MOD2 <-gsub("\\(.*\\)", "", MOD2, perl=TRUE);length(MOD2)
NOUN <-gsub("\\(.*\\)", "", NOUN, perl=TRUE);length(NOUN)
#Get rid of rid of weird punctuation ["()[]\,&}.]
MOD1 <-gsub("[\\\\\\,&}\\.\\(\\)\\\"\\[\\]]", "", MOD1, perl=TRUE)
MOD2 <-gsub("[\\\\\\,&}\\.\\(\\)\\\"\\[\\]]", "", MOD2, perl=TRUE)
NOUN <-gsub("[\\\\\\,&}\\.\\(\\)\\\"\\[\\]]", "", NOUN, perl=TRUE)
#Remove English words and numbers and whitespace from mod1, mod2, noun
MOD1<-gsub("[\\d\\w\\s]", "", MOD1, perl=TRUE)
MOD2<-gsub("[\\d\\w\\s]", "", MOD2, perl=TRUE)
NOUN<-gsub("[\\d\\w\\s]", "", NOUN, perl=TRUE)
#READ BACK IN CODING FILE (END)

# SPLIT INTO MOD AND N ROOTS (START)
#retrieve just the adjective/noun+x roots
MOD1.ROOT<-exact.matches("^[^/]*", MOD1_ANALYSIS, pcre=TRUE, case.sens=FALSE)[[1]];head(MOD1.ROOT)
MOD2.ROOT<-exact.matches("^[^/]*", MOD2_ANALYSIS, pcre=TRUE, case.sens=FALSE)[[1]];head(MOD2.ROOT)
NOUN.ROOT<-exact.matches("^[^/]*", NOUN_ANALYSIS, pcre=TRUE, case.sens=FALSE)[[1]];head(NOUN.ROOT)
#Remove English words and numbers and whitespace from roots
MOD1.ROOT<-gsub("[\\d\\w\\s]", "", MOD1.ROOT, perl=TRUE)
MOD2.ROOT<-gsub("[\\d\\w\\s]", "", MOD2.ROOT, perl=TRUE)
NOUN.ROOT<-gsub("[\\d\\w\\s]", "", NOUN.ROOT, perl=TRUE)
# SPLIT INTO MOD AND N ROOTS (END)

#LENGTH (START)
# compute their lengths
LENGTH1 <- nchar(MOD1) 
LENGTH2 <- nchar(MOD2)
#LENGTH (END)

#SEMOPEN (START)
 #get the freq for paa directly preceding noun, get type counts of those nouns  
 kaist.files <- dir(".", full.names=TRUE) #
 kaist.mod.matches <-c()
 kaist.modn.matches <-c()
for (i in 1:length(kaist.files)){
  cat(i, "\n")
  curr.kaist.file<-scan(kaist.files[i], what=character(0), sep="\n", quiet=TRUE)
  # remove metadata lines
  curr.kaist.file<-gsub("<.*", "", curr.kaist.file, perl=TRUE)
  curr.mod.matches<-exact.matches("([^\\s]+n[^\\s]+하/...(.)?+[^는].etm[^\\s]+xsn(..)?|[^\\s]+paa([^\\s]+etm)?)\\s+", curr.kaist.file, pcre=TRUE, case.sens=FALSE)[[1]]
  curr.modn.matches <-exact.matches("([^\\s]+n[^\\s]+하/...(.)?+[^는].etm|[^\\s]+xsn(..)?|[^\\s]+paa([^\\s]+etm)?)\\s+[^\\s]+\\s+[^\\s]+n[^\\s]*", curr.kaist.file, pcre=TRUE, case.sens=FALSE)[[1]] #this should get adjectives + noun roots
  # add to previous matches of adjective doublets
  kaist.mod.matches <- c(kaist.mod.matches, curr.mod.matches)
  kaist.modn.matches <- c(kaist.modn.matches, curr.modn.matches)
} 
length(kaist.mod.matches);head(kaist.mod.matches)
length(kaist.modn.matches);head(kaist.modn.matches)
#clean up /sl (slang) tags 
kaist.mod.matches <-gsub("(^|\t).?/sl\\+", "", kaist.mod.matches, perl=TRUE)
kaist.modn.matches <-gsub("(^|\t).?/sl\\+", "", kaist.modn.matches, perl=TRUE)
#get rid of anything within parentheses in analysis columns (so roots can be compared)
kaist.mod.matches <-gsub("\\(.*\\)", "", kaist.mod.matches, perl=TRUE);length(kaist.mod.matches)
kaist.modn.matches <-gsub("\\(.*\\)", "", kaist.modn.matches, perl=TRUE);length(kaist.modn.matches)
#Get rid of rid of weird punctuation ["()[]\,&.]
kaist.mod.matches <-gsub("[\\\\\\,&}\\.\\(\\)\\\"\\[\\]]+", "", kaist.mod.matches, perl=TRUE)
kaist.modn.matches <-gsub("[\\\\\\,&}\\.\\(\\)\\\"\\[\\]]+", "", kaist.modn.matches, perl=TRUE)

#for kaist mod.matches:
#get mod root: remove everything past first morpheme
kaist.mod.roots<-gsub("/[^\\s]+", "", kaist.mod.matches, perl=TRUE);head(kaist.mod.roots);length(kaist.mod.roots);length(kaist.mod.matches)
kaist.mod.roots<-gsub("\\t", "", kaist.mod.roots, perl=TRUE);head(kaist.mod.roots)
#clean tags out of kaist.mod.matches
#for kaist.modn.matches: get the adjective and noun roots : 
##split mod and n
kaist.modn.split<-strsplit(kaist.modn.matches, "\\t", perl=TRUE)
##results in three columns, mod analysis, noun, noun analysis
#mod analysis : remove everything past first morpheme
kaist.modn.modroot<-sapply(kaist.modn.split, "[", 1)
kaist.modn.modroot<-gsub("/[^\\s]+", "", kaist.modn.modroot, perl=TRUE);head(kaist.modn.modroot)
kaist.modn.modroot<-gsub("[\\d\\w\\s]", "", kaist.modn.modroot, perl=TRUE)
#noun analysis: remove everything past first morpheme 
kaist.modn.nroot<-sapply(kaist.modn.split, "[", 3) 
kaist.modn.nroot <-gsub("/[^\\s]+", "", kaist.modn.nroot, perl=TRUE);head(kaist.modn.nroot)
kaist.modn.nroot<-gsub("[\\d\\w\\s]", "", kaist.modn.nroot, perl=TRUE)
#paste back together to evaluate uniqueness of combination of two lemmas
kaist.modn.roots<-paste(kaist.modn.modroot, kaist.modn.nroot, sep=" ")
head(kaist.modn.roots);length(kaist.modn.roots)
kaist.modn.unique <-unique(kaist.modn.roots);head(kaist.modn.unique);length(kaist.modn.unique)
#then split again so that roots are only checked for mod
kaist.modn.unique.split<-strsplit(kaist.modn.unique, "\\s", perl=TRUE)
kaist.modn.modroot.unique<-sapply(kaist.modn.unique.split, "[", 1);head(kaist.modn.modroot.unique)
#number of noun types for mod 1 and mod 2
mod1.ntype<-c()
for (j in 1:length(MOD1.ROOT)){ 
		mod1.ntype<-c(mod1.ntype, sum(grepl((MOD1.ROOT[j]), kaist.modn.modroot.unique)))	
	}
head(mod1.ntype);length(mod1.ntype)#3911
hist(mod1.ntype)
plot(sort(table(mod1.ntype)))
#check what is the highest value: 
MOD1.ROOT[which(mod1.ntype==(max(mod1.ntype)))] #같, n=9129 
#MOD 2
mod2.ntype<-c()
for (j in 1:length(MOD2.ROOT)){ 
		mod2.ntype<-c(mod2.ntype, sum(grepl((MOD2.ROOT[j]), kaist.modn.modroot.unique)))	
	}
length(mod2.ntype) #
hist(mod2.ntype)
plot(sort(table(mod2.ntype)))
#check what is the highest value: 
MOD2.ROOT[which(mod2.ntype==(max(mod2.ntype)))] 
#corpus token freq of mods
	##mod 1
mod1.freq<-c()
for (j in 1:length(MOD1.ROOT)){ 
		mod1.freq <-c(mod1.freq, sum(grepl((MOD1.ROOT[j]), kaist.mod.roots)))	
	}
head(mod1.freq);length(mod1.freq)#3911
unique(MOD1.ROOT[which(mod1.freq==(max(mod1.freq)))])
hist(mod1.freq)
plot(sort(table(mod1.freq)))
	##mod 2
mod2.freq<-c()
for (j in 1:length(MOD2.ROOT)){ 
		mod2.freq <-c(mod2.freq, sum(grepl((MOD2.ROOT[j]), kaist.mod.roots)))	
	}
head(mod2.freq);length(mod2.freq)#3911
unique(MOD2.ROOT[which(mod2.freq==(max(mod2.freq)))])
hist(mod2.freq)
plot(sort(table(mod2.freq)))
#Calculate Semantic openness
SEMOPEN1<-mod1.ntype/mod1.freq
SEMOPEN2<-mod2.ntype/mod2.freq
#Prediction SEMOPEN1 < SEMOPEN2 NOUN
SEMOPENDIFF<-SEMOPEN1-SEMOPEN2
##SEM OPEN (END)

#IND COMP (START)
  kaist.files <- dir(".", full.names=TRUE) #
 kaist.comp.matches <-c()
for (i in 1:length(kaist.files)){
  cat(i, "\n")
  curr.kaist.file<-scan(kaist.files[i], what=character(0), sep="\n", quiet=TRUE)
  # remove metadata lines
  curr.kaist.file<-gsub("<.*", "", curr.kaist.file, perl=TRUE) 
  curr.comp.matches<-exact.matches("\\t(보다|제일|더)\\t[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+", curr.kaist.file, pcre=TRUE, case.sens=FALSE)[[1]]
  kaist.comp.matches <- c(kaist.comp.matches, curr.comp.matches)  
}
length(kaist.comp.matches);head(kaist.comp.matches)
kaist.comp.matches2<-exact.matches("(보다|제일|더)\\t[^\\s]+\\s+[^\\s]+\\s+([^\\s]+n[^\\s]+하/...(.)?+[^는].etm|[^\\s]+xsn(..)?|[^\\s]+paa([^\\s]+(etm|ecc))?)", kaist.comp.matches, pcre=TRUE)[[1]]
length(kaist.comp.matches2);head(kaist.comp.matches2)#10,210
#clean up /sl (slang) tags 
kaist.comp.matches2 <-gsub("(^|\t).?/sl\\+", "", kaist.comp.matches2, perl=TRUE)
#get rid of anything within parentheses (so roots can be compared)
kaist.comp.matches2 <-gsub("\\(.*\\)", "", kaist.comp.matches2, perl=TRUE);length(kaist.comp.matches2)
#Get rid of rid of weird punctuation ["()[]]
kaist.comp.matches2 <-gsub("[\\\\\\,&}\\.\\(\\)\\\"\\[\\]]", "", kaist.comp.matches2, perl=TRUE)
length(kaist.comp.matches2);head(kaist.comp.matches2)
#for kaist comp.matches:
#get comp mod root: remove everything past first morpheme
kaist.comp.roots<-exact.matches("([^\\s]+n[^\\s]+하/...(.)?+[^는].etm|[^\\s]+xsn(..)?|[^\\s]+paa([^\\s]+(etm|ecc))?)", kaist.comp.matches2, pcre=TRUE)[[1]];head(kaist.comp.roots);length(kaist.comp.roots)
kaist.comp.roots<-gsub("/[^\\s]+", "", kaist.comp.roots, perl=TRUE);head(kaist.comp.roots);length(kaist.comp.roots)
#Remove English words and numbers and whitespace 
kaist.comp.roots <-gsub("[\\d\\w\\s]", "", kaist.comp.roots, perl=TRUE)
mod1.comp<-c()
for (j in 1:length(MOD1.ROOT)){ 
		mod1.comp<-c(mod1.comp, sum(grepl((MOD1.ROOT[j]), kaist.comp.roots)))	
	}
head(mod1.comp);length(mod1.comp)#3911
hist(mod1.comp)
plot(sort(table(mod1.comp)))
#check what is the highest value: 
unique(MOD1.ROOT[which(mod1.comp==(max(mod1.comp)))])
mod2.comp<-c()
for (j in 1:length(MOD2.ROOT)){ 
		mod2.comp<-c(mod2.comp, sum(grepl((MOD2.ROOT[j]), kaist.comp.roots)))	
	}
head(mod2.comp);length(mod2.comp)#3911
hist(mod2.comp)
plot(sort(table(mod2.comp)))
#check what is the highest value: 
unique(MOD2.ROOT[which(mod2.comp==(max(mod2.comp)))])
INDCOMP1<-mod1.comp/mod1.freq
INDCOMP2<-mod2.comp/mod2.freq
#IND COMP (END)

#SUBJOBJ (START) 
SUBJOBJ1<-as.numeric(MOD1_SUBJOBJ)
SUBJOBJ2<-as.numeric(MOD2_SUBJOBJ)
SUBJOBJ1[which(MOD1=="새로운")]<-c(7)
SUBJOBJ2[which(MOD2=="새로운")]<-c(7)
#SUBJOBJ (END)

#NSPEC FREQ (START)
#f(Adj+N)/f(N)
#Get Noun frequencies from file
kaist.n.freqs<-read.table(file=("../perl/kaist.roots.txt"), quote="", sep="\t", comment.char="")
n.match<-c()
for (j in 1:length(NOUN.ROOT)){ 
		if (sum(kaist.n.freqs$V1==NOUN.ROOT[j])==0){
		n.match<-c(n.match, 0)
		} else {
			n.match <-c(n.match, kaist.n.freqs$V2[which(kaist.n.freqs$V1==NOUN.ROOT[j])])
			}
	}
length(n.match); head(n.match)
NOUN.ROOT[which(n.match==0)]
NOUN_ANALYSIS[which(n.match==0)]
which(n.match==0)
#fix the mis-coded ones and replace derived nouns with their unanalyzed form
#commented out ones were found to be fine by themselves so will search by themselves
NOUN.ROOT[4]<-c("일") #일에대해 
NOUN.ROOT[6]<-c("범위")#범위띠열판 
NOUN.ROOT[13]<-c("에토스")#에토스마음
NOUN.ROOT[89]<-c("저서")#저서우주찬미
NOUN.ROOT[218]<-c("흐름")#흐르
NOUN.ROOT[251]<-c("기쁨") #기쁘, 기쁘/paa+ㅁ/etn+을/jco
NOUN.ROOT[322]<-c("생냉이")#생냉이가
#NOUN.ROOT[343]<-c("")#수연화노리개 "accessory"
#NOUN.ROOT[352]<-c("")#횡혈식 "something"
#NOUN.ROOT[378]<-c("") #박견 "something"
NOUN.ROOT[382]<-c("떨림") #떨리, 떨리/pvg+ㅁ/etn+과/jct
# this needs to be removed NOUN.ROOT[458]<-c("") #따라예까  orig = 따라\예까지
NOUN.ROOT[461]<-c("채찍") #채찍이수없 , orig = 채찍이\수없 , 채찍 'whip'
NOUN.ROOT[462]<-c("하늘") #하늘을내보 , orig = 하늘을\내보인다 , 하늘 'sky, heavens'
NOUN.ROOT[464]<-c("여름") #여름이었다, 여름 'summer'
NOUN.ROOT[470]<-c("그림자") #그림자푸들푸들, orig = 그림자\푸들푸들, 그림자 'shadow'
NOUN.ROOT[475]<-c("것") #것이냐털복숭 , orig = 것이냐\털복숭이 
NOUN.ROOT[476]<-c("사람") #사람들은죽, orig = 사람들은\죽이 
NOUN.ROOT[517]<-c("겨울나기") #겨울나, 겨울나/pvg+기/etn+에/jca
NOUN.ROOT[535]<-c("제품") #제품이다． 제품 'product, goods'
NOUN.ROOT[539]<-c("눈썹") #눈썹，
NOUN.ROOT[663]<-c("밀리") #밀리, 65밀리/ncn
# NOUN.ROOT[668]<-c("") #이륜  this one seems fine, weird it wasn't found before
NOUN.ROOT[703]<-c("젊이") #젊이들 , 젊이 is usually 젊은이 , could search for that separately too
#Run frequency table match again 
missing.n.posit<-c(which(n.match==0))
for (j in missing.n.posit){ 
		if (sum(kaist.n.freqs$V1==NOUN.ROOT[j])==0){
		n.match[j]<-c(0)
		} else {
			n.match[j] <-c(1+kaist.n.freqs$V2[which(kaist.n.freqs$V1==NOUN.ROOT[j])])
			}
	}
length(n.match); head(n.match) #176 1475  153    0   13    0
NOUN.ROOT[which(n.match==0)] 
NOUN_ANALYSIS[which(n.match==0)]
which(n.match==0)
#
kaist.files <- dir(".", full.names=TRUE) #
missing.n.posit<-c(which(n.match==0)) 
missing.n<-c(NOUN.ROOT[which(n.match==0)])
missing.n.freqs<-c(rep(0, length(which(n.match==0))))
missing.n.posit;missing.n;missing.n.freqs
for (i in 1:length(kaist.files)){
  cat(i, "\n")
  curr.kaist.file<-scan(kaist.files[i], what=character(0), sep="\n", quiet=TRUE)
  # remove metadata lines
  curr.kaist.file<-gsub("<.*", "", curr.kaist.file, perl=TRUE)
  for (j in 1:length(missing.n)) {
  	   missing.n.freqs[j] <-c(missing.n.freqs[j]+length(exact.matches(missing.n[j], curr.kaist.file, pcre=TRUE, case.sens=FALSE)[[1]])) 	   
  }
}
missing.n.posit;missing.n;missing.n.freqs
#now replace the n.match values with the new ones
n.match[missing.n.posit]<-c(missing.n.freqs)
head(n.match);length(n.match)#
hist(n.match)
plot(sort(table(n.match)))
max(n.match) #337,070
unique(NOUN.ROOT[which(n.match ==(max(n.match)))]) #것
#p(Adj|N) = p(Adj+N)/p(N) = f(Adj+N)/f(N)

#mod root noun root corpus tokens: kaist.modn.roots
MOD1NROOTS<-paste(MOD1.ROOT, NOUN.ROOT, sep=" ")
MOD2NROOTS<-paste(MOD2.ROOT, NOUN.ROOT, sep=" ")
mod1nmatches <-c()
for (j in 1:length(MOD1NROOTS)){ 
		mod1nmatches <-c(mod1nmatches, sum(grepl(MOD1NROOTS[j], kaist.modn.roots)))	
	}
head(mod1nmatches);length(mod1nmatches)#3911
unique(MOD1NROOTS[which(mod1nmatches ==(max(mod1nmatches)))])
hist(mod1nmatches)
plot(sort(table(mod1nmatches)))
mod2nmatches <-c()
for (j in 1:length(MOD2NROOTS)){ 
		mod2nmatches <-c(mod2nmatches, sum(grepl(MOD2NROOTS[j], kaist.modn.roots)))	
	}
head(mod2nmatches);length(mod2nmatches)#3911
unique(MOD2NROOTS[which(mod2nmatches ==(max(mod2nmatches)))])
hist(mod2nmatches)
plot(sort(table(mod2nmatches)))
NSPEC1<-(mod1nmatches/n.match) #one word, 아프, has value higher than 1
NSPEC2<-(mod2nmatches/n.match)
#NSPEC FREQ (END)

#GENFREQ (START)
FREQ1<-mod1.freq
FREQ2<-mod2.freq
#GENFREQ (END)

#NOMCHAR (START)
#English : freq of adj tag/ (freq of adj tag+ freq of noun tag)
#more noun like (lower value of NOMCHAR) goes closer to noun
#Prediction:NOMCHAR 1 > NOMCHAR 2 NOUN 
#n in root vs. paa in root
#There are no 'xsn' matches in MOD1 and only 8 in MOD2 so contrasting just that w/stative verbs seems no good
#optionally modifiers with paa in root but with 하  could be a middle level of nouniness
#so instead find first morpheme tag
sum(grepl("^\\W+/n", MOD1_ANALYSIS, perl=TRUE)) #0
sum(grepl("^\\W+/p", MOD1_ANALYSIS, perl=TRUE)) #711
sum(grepl("^\\W+/n", MOD2_ANALYSIS, perl=TRUE)) #8
sum(grepl("^\\W+/p", MOD2_ANALYSIS, perl=TRUE)) #703
NOMCHAR1<-grepl("^\\W+/n", MOD1_ANALYSIS, perl=TRUE)
NOMCHAR2<-grepl("^\\W+/n", MOD2_ANALYSIS, perl=TRUE)
NOMCHAR1[which(NOMCHAR1=="FALSE")]<-c("adj")
NOMCHAR1[which(NOMCHAR1=="TRUE")]<-c("n")
NOMCHAR2[which(NOMCHAR2=="FALSE")]<-c("adj")
NOMCHAR2[which(NOMCHAR2=="TRUE")]<-c("n")
table(NOMCHAR1);table(NOMCHAR2)
#NOMCHAR (END)

