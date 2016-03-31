#rm(list=ls(all=TRUE))

#Loading script for nice regex searching with greedy return of string matches
source("http://www.linguistics.ucsb.edu/faculty/stgries/exact_matches.r")

# #Check # of tokens in the corpus
 # kaist.files <- dir(".", full.names=TRUE) #
 # kaist.word.counts <-c()
# for (i in 1:length(kaist.files)){
  # cat(i, "\n")
  # curr.kaist.file<-scan(kaist.files[i], what=character(0), sep="\n", quiet=TRUE)
  # # remove metadata lines
  # curr.kaist.file<-gsub("<.*", "", curr.kaist.file, perl=TRUE)
   # curr.word.count <-length(exact.matches("[^\\s]+\\s+[^\\s]+\\/[^\\s]+", curr.kaist.file, pcre=TRUE, case.sens=FALSE)[[1]])
  # # add to previous matches of adjective doublets 
  # kaist.word.counts <- c(kaist.word.counts, curr.word.count)
# } 
filewc<-read.table(file="../wc.txt", header=FALSE, sep="", quote="")
str(filewc)
filewc<-c(filewc$V2)
filewc<-as.numeric(filewc)
sum(filewc) #token count = 61,060,421
#Their website says the raw corpus is 70 million "phrases" (http://semanticweb.kaist.ac.kr/home/index.php/KAIST_Corpus). However, the token count is really about half that amount, due to wc including the analysis and word form separately. This can be shown by looking at the first, small, file
biblefile<-scan(file="./BIBLE-1.txt.new.pos2", what=character(0), sep="\n", quiet=TRUE) #wc = 272
length(exact.matches("[^\\s]+\\s+[^\\s]+\\/[^\\s]+", biblefile, pcre=TRUE, case.sens=FALSE)[[1]]) #111
length(exact.matches("\\t", biblefile, pcre=TRUE, case.sens=FALSE)[[1]]) #212
#so really it's about 30 million words 

kaist.analyzed <-read.table(file=("../KAIST_raw_pos_modmodn.csv"), quote="", sep="\t", header=TRUE, comment.char="")
# str(kaist.analyzed)
kaist.coded <-read.table(file=("../KAIST_Wona_Coding_final.csv"), quote="", sep="\t", header=TRUE, comment.char="")
# str(kaist.coded)
# #attach both bc the coding columns are not in kaist.analyzed so won't overwrite with empty
attach(kaist.coded)
attach(kaist.analyzed)
#add analysis columns back in from KAIST raw mod mod n
#kaist.mod.mod.n<-data.frame(CASE, FILE, GOOD, MOD1, MOD1_ANALYSIS, MOD1_DEF,MOD1_SUBJOBJ, MOD1_AFFLOAD, MOD2, MOD2_ANALYSIS, MOD2_DEF, MOD2_SUBJOBJ, MOD2_AFFLOAD, NOUN, NOUN_ANALYSIS,NOUN_DEF)
# #Stuff to make the Definition variable is below
# GOOD<-kaist.coded$GOOD
# Definition <-kaist.coded$Definition
# M1_SO<-kaist.coded$MOD1_SUBJOBJ
# M2_SO<-kaist.coded$MOD2_SUBJOBJ
# M1_AL<-kaist.coded$MOD1_AFFLOAD
# M2_AL<-kaist.coded$MOD2_AFFLOAD
# kaist.mod.mod.n<-data.frame(GOOD, Definition, M1_SO, M2_SO, M1_AL, M2_AL)
# # str(kaist.mod.mod.n)
# detach(kaist.coded)
# detach(kaist.analyzed)
# good.matchn<-grepl("n", kaist.mod.mod.n$GOOD, perl=TRUE, ignore.case=TRUE)
# kaist.mod.mod.n<-kaist.mod.mod.n[which(good.matchn=="FALSE"),]
# library(gdata)
# kaist.mod.mod.n<-drop.levels(kaist.mod.mod.n)
# good.matchcheck<-grepl("check", kaist.mod.mod.n$GOOD, perl=TRUE, ignore.case=TRUE)
# kaist.mod.mod.n<-kaist.mod.mod.n[which(good.matchcheck=="FALSE"),]
# kaist.mod.mod.n<-drop.levels(kaist.mod.mod.n)
# subjobj1.match<-grepl("\\d", kaist.mod.mod.n$M1_SO, perl=TRUE, ignore.case=TRUE) 
# kaist.mod.mod.n<-kaist.mod.mod.n[which(subjobj1.match=="TRUE"),]
# kaist.mod.mod.n<-drop.levels(kaist.mod.mod.n)
# subjobj2.match<-grepl("\\d", kaist.mod.mod.n$M2_SO, perl=TRUE, ignore.case=TRUE)
# kaist.mod.mod.n<-kaist.mod.mod.n[which(subjobj2.match=="TRUE"),]
# kaist.mod.mod.n<-drop.levels(kaist.mod.mod.n)
# affload1.match<-grepl("\\d", kaist.mod.mod.n$M1_AL, perl=TRUE, ignore.case=TRUE)
# kaist.mod.mod.n<-kaist.mod.mod.n[which(affload1.match=="TRUE"),]
# kaist.mod.mod.n<-drop.levels(kaist.mod.mod.n)
# affload2.match<-grepl("\\d", kaist.mod.mod.n$M2_AL, perl=TRUE, ignore.case=TRUE)
# kaist.mod.mod.n<-kaist.mod.mod.n[which(affload2.match=="TRUE"),]
# kaist.mod.mod.n2<-drop.levels(kaist.mod.mod.n)
# str(kaist.mod.mod.n2) #1810 obs 
# Definition<-kaist.mod.mod.n2$Definition

# attach(kaist.mod.mod.n)
# sum(grepl("\\(|\\)", MOD1, perl=TRUE))
# sum(grepl("\\(|\\)", MOD2, perl=TRUE))
# sum(grepl("\\(|\\)", NOUN, perl=TRUE))
# #kaist.mod.mod.n<-kaist.mod.mod.n[grepl(kaist.mod.mod.n$GOOD)]
# #clean up /sl (slang) tags 
# MOD1_ANALYSIS <-gsub("(^|\t).?/sl\\+", "", MOD1_ANALYSIS, perl=TRUE)
# MOD2_ANALYSIS <-gsub("(^|\t).?/sl\\+", "", MOD2_ANALYSIS, perl=TRUE)
# NOUN_ANALYSIS <-gsub("(^|\t).?/sl\\+", "", NOUN_ANALYSIS, perl=TRUE)
# #get rid of anything within parentheses in analysis columns (so roots can be compared)
# MOD1_ANALYSIS <-gsub("\\(.*\\)", "", MOD1_ANALYSIS, perl=TRUE);length(MOD1_ANALYSIS)
# MOD2_ANALYSIS <-gsub("\\(.*\\)", "", MOD2_ANALYSIS, perl=TRUE);length(MOD2_ANALYSIS)
# NOUN_ANALYSIS <-gsub("\\(.*\\)", "", NOUN_ANALYSIS, perl=TRUE);length(NOUN_ANALYSIS)
# #or on one side of parentheses
# MOD1_ANALYSIS <-gsub("\\([^\\)]+$", "", MOD1_ANALYSIS, perl=TRUE);length(MOD1_ANALYSIS)
# MOD2_ANALYSIS <-gsub("\\([^\\)]+$", "", MOD2_ANALYSIS, perl=TRUE);length(MOD2_ANALYSIS)
# NOUN_ANALYSIS <-gsub("\\([^\\)]+$", "", NOUN_ANALYSIS, perl=TRUE);length(NOUN_ANALYSIS)
# #Get rid of rid of weird punctuation ["()[]\}&,]
# MOD1_ANALYSIS <-gsub("[\\\\\\,&}\\.\\(\\)\\\"\\[\\]]", "", MOD1_ANALYSIS, perl=TRUE)
# MOD2_ANALYSIS <-gsub("[\\\\\\,&}\\.\\(\\)\\\"\\[\\]]", "", MOD2_ANALYSIS, perl=TRUE)
# NOUN_ANALYSIS <-gsub("[\\\\\\,&}\\.\\(\\)\\\"\\[\\]]", "", NOUN_ANALYSIS, perl=TRUE)
# #clean up /sl (slang) tags 
# MOD1 <-gsub("(^|\t).?/sl\\+", "", MOD1, perl=TRUE)
# MOD2 <-gsub("(^|\t).?/sl\\+", "", MOD2, perl=TRUE)
# NOUN <-gsub("(^|\t).?/sl\\+", "", NOUN, perl=TRUE)
# #get rid of anything within parentheses in analysis columns (so roots can be compared)
# MOD1 <-gsub("\\(.*\\)", "", MOD1, perl=TRUE);length(MOD1)
# MOD2 <-gsub("\\(.*\\)", "", MOD2, perl=TRUE);length(MOD2)
# NOUN <-gsub("\\(.*\\)", "", NOUN, perl=TRUE);length(NOUN)
# #Get rid of rid of weird punctuation ["()[]\,&}.]
# MOD1 <-gsub("[\\\\\\,&}\\.\\(\\)\\\"\\[\\]]", "", MOD1, perl=TRUE)
# MOD2 <-gsub("[\\\\\\,&}\\.\\(\\)\\\"\\[\\]]", "", MOD2, perl=TRUE)
# NOUN <-gsub("[\\\\\\,&}\\.\\(\\)\\\"\\[\\]]", "", NOUN, perl=TRUE)
# #Remove English words and numbers and whitespace from mod1, mod2, noun
# MOD1<-gsub("[\\d\\w\\s]", "", MOD1, perl=TRUE)
# MOD2<-gsub("[\\d\\w\\s]", "", MOD2, perl=TRUE)
# NOUN<-gsub("[\\d\\w\\s]", "", NOUN, perl=TRUE)
# #READ BACK IN CODING FILE (END)

# # SPLIT INTO MOD AND N ROOTS (START)
# #retrieve just the adjective/noun+x roots
# MOD1.ROOT<-exact.matches("^[^/]*", MOD1_ANALYSIS, pcre=TRUE, case.sens=FALSE)[[1]];head(MOD1.ROOT)
# MOD2.ROOT<-exact.matches("^[^/]*", MOD2_ANALYSIS, pcre=TRUE, case.sens=FALSE)[[1]];head(MOD2.ROOT)
# NOUN.ROOT<-exact.matches("^[^/]*", NOUN_ANALYSIS, pcre=TRUE, case.sens=FALSE)[[1]];head(NOUN.ROOT)
# #Remove English words and numbers and whitespace from roots
# MOD1.ROOT<-gsub("[\\d\\w\\s]", "", MOD1.ROOT, perl=TRUE)
# MOD2.ROOT<-gsub("[\\d\\w\\s]", "", MOD2.ROOT, perl=TRUE)
# NOUN.ROOT<-gsub("[\\d\\w\\s]", "", NOUN.ROOT, perl=TRUE)
# # # SPLIT INTO MOD AND N ROOTS (END)
# kaist.mod.mod.n<-data.frame(GOOD, FILE, MOD1, MOD1.ROOT, MOD1_ANALYSIS, MOD1_DEF,MOD1_SUBJOBJ, MOD1_AFFLOAD, MOD2, MOD2.ROOT, MOD2_ANALYSIS, MOD2_DEF, MOD2_SUBJOBJ, MOD2_AFFLOAD, NOUN, NOUN.ROOT, NOUN_ANALYSIS,NOUN_DEF)
# write.table(kaist.mod.mod.n, file="../KAIST.mod.mod.n.csv", quote=F, sep="\t")
kaist.mod.mod.n<-read.table(file=("../KAIST.mod.mod.n.csv"), quote="", sep="\t", header=TRUE, comment.char="")
attach(kaist.mod.mod.n)

#LENGTH (START)
# compute their lengths
RT.LENGTH1 <- nchar(as.vector(MOD1.ROOT))
RT.LENGTH2 <- nchar(as.vector(MOD2.ROOT))
MOD2[949];MOD2.ROOT[949]
RT.LENGTH2[949] <-4#fix this 차이나식 
#LENGTH (END)

#SEMOPEN (START)
#get freq for modifiers directly preceding noun, get type counts of those nouns  
#get the set of adjectives in attributive constructions (adj + n) from corpus
 kaist.files <- dir(".", full.names=TRUE) #
 kaist.mod.matches <-c()
 kaist.modn.matches <-c()
for (i in 1:length(kaist.files)){
  cat(i, "\n")
  curr.kaist.file<-scan(kaist.files[i], what=character(0), sep="\n", quiet=TRUE)
  # remove metadata lines\\\\\\\\\\\\\\\
   curr.mod.matches<-exact.matches("[^\\s]+paa([^\\s]+etm)?", curr.kaist.file, pcre=TRUE, case.sens=FALSE)[[1]] #this has 345292 which is less than modn matches, doesn't make sense, should just get modn matches then extract mod from it
  curr.modn.matches <-exact.matches("[^\\s]+paa([^\\s]+etm)?\\s+[^\\s]+\\s+[^\\s]+n[^\\s]*", curr.kaist.file, pcre=TRUE, case.sens=FALSE)[[1]] #this should get adjectives + noun roots
  # add to previous matches of adjective doublets 
  kaist.mod.matches <- c(kaist.mod.matches, curr.mod.matches)
  kaist.modn.matches <- c(kaist.modn.matches, curr.modn.matches)
} 
length(kaist.mod.matches);head(kaist.mod.matches)
length(kaist.modn.matches);head(kaist.modn.matches)#348142
#write.table(kaist.modn.matches, file="../KAIST_modnmatches_orig.csv", quote=F, sep="\t", row.names=F, col.names=F)
kaist.modn.matches <-read.table(file=("../KAIST_modnmatches_orig.csv"), quote="", sep="\t", header=F, comment.char="")

#need to fix the clean up code below if running on KAIST_modnmatches_orig.csv, because now there are separate columns to deal with
kaist.modn.modroot<-kaist.modn.matches$V1
kaist.modn.nroot<-kaist.modn.matches$V3 # $V2 is the noun full form , $V3 is the noun analysis
kaist.modn.modroot <-gsub(".\\/sl\\+", "", kaist.modn.modroot, perl=TRUE) #clean up /sl (slang) tags
kaist.modn.modroot <-gsub("\\+.\\/sr", "", kaist.modn.modroot, perl=TRUE) #clean up /sl (slang) tags
kaist.modn.modroot <-gsub("\\(.*\\)", "", kaist.modn.modroot, perl=TRUE);length(kaist.modn.modroot) #get rid of anything within parentheses in analysis columns
kaist.modn.nroot <-gsub("[\\\\\\,&}\\.\\(\\)\\\"\\[\\]]+", "", kaist.modn.nroot, perl=TRUE) #Get rid of rid of weird punctuation ["()[]\,&.]
kaist.modn.nroot <-gsub(".\\/sl\\+", "", kaist.modn.nroot, perl=TRUE) #clean up /sl (slang) tags
kaist.modn.nroot <-gsub("\\+.\\/sr", "", kaist.modn.nroot, perl=TRUE) #clean up /sl (slang) tags
kaist.modn.nroot <-gsub("\\(.*\\)", "", kaist.modn.nroot, perl=TRUE);length(kaist.modn.nroot) #get rid of anything within parentheses in analysis columns
kaist.modn.nroot <-gsub("[\\\\\\,&}\\.\\(\\)\\\"\\[\\]]+", "", kaist.modn.nroot, perl=TRUE) #Get rid of rid of weird punctuation ["()[]\,&.]

#get the adjective and noun roots : 
##split mod and n
# kaist.modn.onetab.matches<-c(exact.matches("^[^\t]+\t[^\t]+$",kaist.modn.matches, pcre=TRUE)[[1]])
# kaist.modn.twotab.matches<-c(exact.matches(".+\t.+\t.+",kaist.modn.matches, pcre=TRUE)[[1]])
# kaist.modn.onetab.split<-strsplit(kaist.modn.onetab.matches, "\\t", perl=TRUE)
# kaist.modn.twotab.split<-strsplit(kaist.modn.twotab.matches, "\\t", perl=TRUE)

# ##results in three columns, mod analysis, noun, noun analysis
# #mod analysis : remove everything past first morpheme
# kaist.modn.modroot<-c(sapply(kaist.modn.onetab.split, "[", 1), sapply(kaist.modn.twotab.split, "[", 1))
kaist.modn.modroot<-gsub("\\/[^\\s]+", "", kaist.modn.modroot, perl=TRUE);head(kaist.modn.modroot)
kaist.modn.modroot<-gsub("[\\d\\w\\s]", "", kaist.modn.modroot, perl=TRUE)
#noun analysis: remove everything past first morpheme 
# kaist.modn.nroot2<-sapply(kaist.modn.twotab.split, "[", 3)
# kaist.modn.nroot1<-sapply(kaist.modn.onetab.split, "[", 2)
# kaist.modn.nroot<-c(kaist.modn.nroot1, kaist.modn.nroot2)
#make special set for ncn ncn compounds
# kaist.modn.modncn<-kaist.modn.modroot[grep("ncn\\+[^들이]+\\/ncn", kaist.modn.nroot, perl=TRUE)]
kaist.modn.ncn.posit<-grep("ncn\\+[^들이]+\\/ncn", kaist.modn.nroot, perl=TRUE)  #get positions of compound nouns 
bothncn<-exact.matches("^.+ncn\\+[^들이]+\\/ncn", as.vector(kaist.modn.nroot), pcre=TRUE, case.sens=FALSE)[[1]];length(bothncn) #extract noun compounds
kaist.modn.nroot <-gsub("\\/.+", "", kaist.modn.nroot, perl=TRUE);head(kaist.modn.nroot) #remove everything past first morpheme from full list of noun roots
kaist.modn.nroot<-gsub("[\\d\\w\\s]", "", kaist.modn.nroot, perl=TRUE)
ncn1<-exact.matches("^[^\\/]+", bothncn, pcre=TRUE, case.sens=FALSE)[[1]];head(ncn1);length(ncn1)
ncn2<-exact.matches("ncn\\+[^들이]+\\/ncn", bothncn, pcre=TRUE, case.sens=FALSE)[[1]];head(ncn2);length(ncn2)
ncn2<-gsub("[\\d\\w\\s\\/\\+]", "", ncn2, perl=TRUE);head(ncn2);length(ncn2)
ncnlemmas<-paste(ncn1, ncn2, sep="");head(ncnlemmas);length(ncnlemmas)
kaist.modn.nroot[kaist.modn.ncn.posit]<-ncnlemmas #replace kaist.modn.nroots that are really compound nouns with the original ncn matches
#paste back together to evaluate uniqueness of combination of two lemmas
kaist.modn.roots<-paste(kaist.modn.modroot, kaist.modn.nroot, sep=" ")
head(kaist.modn.roots);length(kaist.modn.roots)
kaist.modn.unique<-unique(kaist.modn.roots);head(kaist.modn.unique);length(kaist.modn.unique) #117202
#then split again so that we can just count the number of times the adj occurs as the number of noun types
kaist.modn.unique.split<-strsplit(kaist.modn.unique, "\\s", perl=TRUE)
kaist.modn.modroot.unique<-sapply(kaist.modn.unique.split, "[", 1);head(kaist.modn.modroot.unique)

#NTYPE for MOD1
mod1.ntype<-c()
for (j in 1:length(MOD1.ROOT)){ 
		mod1.ntype<-c(mod1.ntype, length(which(kaist.modn.modroot.unique==MOD1.ROOT[j])))	
	}
head(mod1.ntype);length(mod1.ntype)
hist(mod1.ntype)
plot(sort(table(mod1.ntype)))

#Add the counts from 1st mod in modmod position
# modmod1n.unique<-unique(paste(MOD1.ROOT, NOUN.ROOT, sep=" "))
# modmod1n.unique.split<-strsplit(modmod1n.unique, "\\s", perl=TRUE)
# modmod1n.modroot.unique<-sapply(modmod1n.unique.split, "[", 1);head(modmod1n.modroot.unique);length(modmod1n.modroot.unique)
# mod1.ntypefrommodmod<-c()
# for (j in 1:length(MOD1.ROOT)){ 
		# mod1.ntypefrommodmod<-c(mod1.ntypefrommodmod, length(which(modmod1n.modroot.unique==MOD1.ROOT[j])))	
	# }
# length(mod1.ntypefrommodmod);head(mod1.ntypefrommodmod)
# mod1.ntype<-c(mod1.ntype+mod1.ntypefrommodmod)
# head(mod1.ntype);length(mod1.ntype)#1810
# write.table(mod1.ntype, file="../KAIST_mod1.ntype.csv", quote=F, row.names=F, col.names=F, sep="\t")
mod1.ntype<-read.table(file="../KAIST_mod1.ntype.csv", quote="", header=F, sep="\t")
mod1.ntype<-mod1.ntype$V1
# which(mod1.ntype==0)# 329  918 1053

#NTYPE for MOD2
# mod2.ntype<-c()
# for (j in 1:length(MOD2.ROOT)){ 
		# mod2.ntype<-c(mod2.ntype, length(which(kaist.modn.modroot.unique==MOD2.ROOT[j])))	
	# }
# length(mod2.ntype) #
# hist(mod2.ntype)
# plot(sort(table(mod2.ntype)))
# #Add the counts from 2nd mod in 1st mod mod position
# mod2.ntypefrommodmod<-c()
# for (j in 1:length(MOD2.ROOT)){ 
		# mod2.ntypefrommodmod<-c(mod2.ntypefrommodmod, length(which(modmod1n.modroot.unique==MOD2.ROOT[j])))	
	# }
# length(mod2.ntypefrommodmod);head(mod2.ntypefrommodmod)
# mod2.ntype<-c(mod2.ntype+mod2.ntypefrommodmod)
# head(mod2.ntype);length(mod2.ntype)#3922
# write.table(mod2.ntype, file="../KAIST_mod2.ntype.csv", quote=F, row.names=F, col.names=F, sep="\t")
mod2.ntype<-read.table(file="../KAIST_mod2.ntype.csv", quote="", header=F, sep="\t")
mod2.ntype<-mod2.ntype$V1
# #check what is the highest value: 
# MOD2.ROOT[which(mod2.ntype==(max(mod2.ntype)))] 

#MOD token frequency
	# # MOD2 Frequency
mod1.freq<-c()
for (j in 1:length(MOD1.ROOT)){ 
		mod1.freq <-c(mod1.freq, length(which(kaist.modn.modroot==MOD1.ROOT[j])))	
	}
#plus 1st position mod in mod mod n
head(mod1.freq);length(mod1.freq)#1810
unique(MOD1.ROOT[which(mod1.freq==(max(mod1.freq)))])
hist(mod1.freq)
plot(sort(table(mod1.freq)))
mod1.freq2<-c()
for (j in 1:length(MOD1.ROOT)){ 
		mod1.freq2 <-c(mod1.freq2, length(which(MOD1.ROOT==MOD1.ROOT[j])))	
	}
head(mod1.freq2);length(mod1.freq2)#1810
unique(MOD1.ROOT[which(mod1.freq2==(max(mod1.freq2)))])
mod1.freq<-c(mod1.freq+mod1.freq2)
#write.table(mod1.freq, file="../KAIST_mod1.freq.csv", row.names=F, col.names=F,quote=F, sep="\t")
mod1.freq<-read.table(file="../KAIST_mod1.freq.csv", header=FALSE, sep="\t", quote="")
mod1.freq<-c(mod1.freq$V1)
	# # MOD2 Frequency
# mod2.freq<-c()
# for (j in 1:length(MOD2.ROOT)){ 
		# mod2.freq <-c(mod2.freq, length(which(kaist.modn.modroot==MOD2.ROOT[j])))	
	# }
# head(mod2.freq);length(mod2.freq)#3911
# unique(MOD2.ROOT[which(mod2.freq==(max(mod2.freq)))])
# hist(mod2.freq)
# plot(sort(table(mod2.freq)))
# mod1root<-as.vector(MOD1.ROOT)
# mod2root<-as.vector(MOD2.ROOT)
# mod2.freq2<-c()
# for (j in 1:length(mod2root)){ 
		# mod2.freq2 <-c(mod2.freq2, length(which(mod1root==mod2root[j])))	
		# }
# head(mod2.freq2);length(mod2.freq2)#2820
# unique(MOD2.ROOT[which(mod2.freq2==(max(mod2.freq2)))])
# mod2.freq<-c(mod2.freq+mod2.freq2)
# mod2.freq[949]<-4 #fix this 차이나식 one
# write.table(mod2.freq, file="../KAIST_mod2.freq.csv", quote=F, row.names=F, col.names=F, sep="\t")
mod2.freq<-read.table(file="../KAIST_mod2.freq.csv", header=FALSE, sep="\t", quote="")
mod2.freq<-c(mod2.freq$V1)
#Calculate Semantic openness
SEMOPEN1<-c(mod1.ntype/mod1.freq)
SEMOPEN2<-c(mod2.ntype/mod2.freq)
# #Prediction SEMOPEN1 > SEMOPEN2 NOUN
SEMOPENDIFF<-c(SEMOPEN1-SEMOPEN2)
#semopen<-data.frame(MOD1, SEMOPEN1, MOD2, SEMOPEN2, SEMOPENDIFF, NOUN) 
#write.table(semopen, file="../KAIST_semopen.csv", quote=F, sep="\t", row.names=F)
#semopen<-read.table(file="../KAIST_semopen.csv", header=T, quote="", sep="\t")
##SEM OPEN (END)

#IND COMP (START)
 kaist.files <- dir(".", full.names=TRUE) #
 kaist.comp.matches <-c()
for (i in 1:length(kaist.files)){
  cat(i, "\n")
  curr.kaist.file<-scan(kaist.files[i], what=character(0), sep="\n", quiet=TRUE)
  # remove metadata lines
  curr.kaist.file<-gsub("<.*", "", curr.kaist.file, perl=TRUE) 
  curr.comp.matches<-exact.matches("([^\\s]+\\s+[^\\s]+\\s+)?[^\\s]+\\s+[^\\s]+\\s+(보다|제일|더)\\t[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+\\s+[^\\s]+", curr.kaist.file, pcre=TRUE, case.sens=FALSE)[[1]]
  kaist.comp.matches <- c(kaist.comp.matches, curr.comp.matches)  
}
length(kaist.comp.matches);head(kaist.comp.matches)#33980
write.table(kaist.comp.matches, file="../KAIST_compmatches.csv", quote=F, sep="\t", row.names=F, col.names=F)
kaist.comp.matches2<-exact.matches("(보다|제일|더)\\t[^\\s]+\\s+[^\\s]+\\s+([^\\s]+xsn(..)?|[^\\s]+paa[^\\s]+|[^\\s]+\\s+[^\\s]+\\s+([^\\s]+xsn(..)?|[^\\s]+paa[^\\s]+))", kaist.comp.matches, pcre=TRUE)[[1]] # ONLY paa and xsn because that's all we have in kaistresults, and also don't care what is after paa because it seems to be consistent that what is after (보다|제일|더) must be the "more than" thing , but also can allow one non paa word in between, 것보다\t것/nbn+보다/jca\t훨씬\t훨씬/mag\t더\t더/mag\t하나님께서\t하나님께/ncn+서/jca\t기뻐하시는\t기쁘/paa+어/ecx+하/px+시/ep+는/etm"
length(kaist.comp.matches2);head(kaist.comp.matches2)#11134
#clean up /sl (slang) tags 
kaist.comp.matches2 <-gsub("(^|\t).?/sl\\+", "", kaist.comp.matches2, perl=TRUE)
#get rid of anything within parentheses (so roots can be compared)
kaist.comp.matches2 <-gsub("\\(.*\\)", "", kaist.comp.matches2, perl=TRUE);length(kaist.comp.matches2)
#Get rid of rid of weird punctuation ["()[]]
kaist.comp.matches2 <-gsub("[\\\\\\,&}\\.\\(\\)\\\"\\[\\]]", "", kaist.comp.matches2, perl=TRUE)
length(kaist.comp.matches2);head(kaist.comp.matches2)
write.table(kaist.comp.matches2, file="../KAIST_compmatches2.csv", quote=F, sep="\t", row.names=F, col.names=F)
#for kaist comp.matches:
#get comp mod root: get the 
kaist.comp.roots<-exact.matches("[^\\s]+xsn(..)?|[^\\s]+paa[^\\s]+", kaist.comp.matches2, pcre=TRUE)[[1]];head(kaist.comp.roots);length(kaist.comp.roots)
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
indcomp<-data.frame(MOD1, INDCOMP1, MOD2, INDCOMP2, NOUN)
write.table(indcomp, file="../KAIST_indcomp.csv", quote=F, sep="\t")
indcomp<-read.table(file="../KAIST_indcomp.csv", header=T, quote="", sep="\t")
#IND COMP (END)


# #NSPEC FREQ (START)
# #f(Adj+N)/f(N)

##fix NOUN.ROOT for compound nouns
NOUN_ANALYSIS<-as.vector(NOUN_ANALYSIS)
NOUN.ROOT<-as.vector(NOUN.ROOT)
NOUN_ANALYSIS[grep("ncn\\+[^들이]+\\/ncn", NOUN_ANALYSIS, perl=TRUE)]
ncnposit<-c(grep("ncn\\+[^들이]+\\/ncn", NOUN_ANALYSIS, perl=TRUE))
ncn_ncn<-c(exact.matches("[^\\/]+\\/ncn\\+[^들이]+\\/ncn", NOUN_ANALYSIS, pcre=TRUE, case.sens=FALSE)[[1]])
NOUN.ROOT2<-NOUN.ROOT
ncn1<-exact.matches("^[^\\/]+", ncn_ncn, pcre=TRUE, case.sens=FALSE)[[1]];head(ncn1);length(ncn1)
ncn2<-exact.matches("ncn\\+[^들이]+\\/ncn", ncn_ncn, pcre=TRUE, case.sens=FALSE)[[1]];head(ncn2);length(ncn2)
ncn2<-gsub("[\\d\\w\\s\\/\\+]", "", ncn2, perl=TRUE);head(ncn2);length(ncn2)
ncnlemmas<-paste(ncn1, ncn2, sep="");head(ncnlemmas);length(ncnlemmas)
NOUN.ROOT2[ncnposit]<-ncnlemmas
# #Get Noun frequencies from file
# kaist.n.freqs<-read.table(file=("../perl/kaist.roots3.txt"), quote="", sep="\t", comment.char="")
# #match up with NOUN.ROOT
# n.match<-c()
# for (j in 1:length(NOUN.ROOT2)){ 
		# if (sum(kaist.n.freqs$V1==NOUN.ROOT2[j])==0){
		# n.match<-c(n.match, 0)
		# } else {
			# n.match <-c(n.match, kaist.n.freqs$V2[which(kaist.n.freqs$V1==NOUN.ROOT2[j])])
			# }
	# }
# length(n.match); head(n.match)
# NOUN.ROOT2[which(n.match==0)]
# NOUN_ANALYSIS[which(n.match==0)]
# which(n.match==0)
# n.match2<-n.match

# #need to add counts for some mistagged things 
# #n.match2[4]<-c("일") #일에대해 
# n.match2[4]<-c(sum(n.match2[4], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="일")]))
# #NOUN.ROOT[83]<-c("저서")#저서우주찬미 
# n.match2[83]<-c(sum(n.match2[83], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="저서")]));n.match2[83]

# #NOUN.ROOT[211]<-c("흐르")#흐르
	# #need to add results for 흐름 separately 
# n.match2[211];n.match2[211]<-c(sum(n.match2[211], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="흐름")]));n.match2[211]

# #NOUN.ROOT[244]<-c("기쁨")#기쁘, 기쁘/paa+ㅁ/etn+을/jco
	# #need to add results for "기쁨" separately 
# n.match2[244];n.match2[244]<-c(sum(n.match2[244], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="기쁨")]));n.match2[244]

# #NOUN.ROOT[315]<-c("생냉이")#생냉이가,
# n.match2[315];n.match2[315]<-c(sum(n.match2[315], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="생냉이")]));n.match2[315] 

# #NOUN.ROOT[374]<-c("떨림")#떨리/pvg+ㅁ/etn+과/jct" 
# n.match2[374];n.match2[374]<-c(sum(n.match2[374], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="떨림")]));n.match2[374]

# #NOUN.ROOT[536] #"", "/sr+음식/ncn+점/ncn"
# n.match2[536];n.match2[536]<-c(sum(kaist.n.freqs$V2[which(kaist.n.freqs$V1=="음식점")]));n.match2[536] 
# #UPDATE: will be updated and added to kaist.roots_ncn.txt : 음식/ncn+점/ncn 	469

# #NOUN.ROOT[900] #기쁘/paa+ㅁ/etn"
# n.match2[900];n.match2[900]<-c(sum(n.match2[900], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="기쁨")]));n.match2[900]  
   
# #NOUN.ROOT[1037]#꽃들이오무라진/ncn"  add "꽃"
# n.match2[1037];n.match2[1037]<-c(sum(n.match2[1037], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="꽃")]));n.match2[1037] 

# #NOUN.ROOT[1090]#애를…/ncn+/sl"  
# n.match2[1090];n.match2[1090]<-c(sum(n.match2[1090], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="애")]));n.match2[1090] 

# #NOUN.ROOT[1091] #아가씨였는데…/ncn+/sl"  
# n.match2[1091];n.match2[1091]<-c(sum(n.match2[1091], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="아가씨")]));n.match2[1091]

# NOUN.ROOT[1145] #비가\무너진,  noun is "비가", root "비" 'rain'
# n.match2[1145];n.match2[1145]<-c(sum(n.match2[1145], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="비")]));n.match2[1145]

# NOUN.ROOT[1149] #하늘에\새하얀 # 
# n.match2[1149];n.match2[1149]<-c(sum(n.match2[1149], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="하늘")]));n.match2[1149]

# #NOUN.ROOT[1151] #새투명한, orig "새,\투명한" , root 새 'bird'
# n.match2[1151];n.match2[1151]<-c(sum(n.match2[1151], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="새")]));n.match2[1151])
# n.match2[1151]<-c(n.match2[1151]+1);n.match2[1151]# adding the orig form itself

# #NOUN.ROOT[1153] #채찍이수없/ncn+이/jcc"  , orig = 채찍이\수없 , 채찍 'whip'
# n.match2[1153];n.match2[1153]<-c(sum(n.match2[1153], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="채찍")]));n.match2[1153]

# NOUN.ROOT[1155]<-c("하늘") #하늘을내보/ncn+이/jp+ㄴ다/ef", 하늘을내보 , orig = 하늘을\내보인다 , 하늘 'sky, heavens'
# n.match2[1155];n.match2[1155]<-c(sum(n.match2[1155], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="하늘")]));n.match2[1155]

# NOUN.ROOT[1157] #시절\50원/ncn+을/jco" 
# n.match2[1157];n.match2[1157]<-c(sum(n.match2[1157], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="시절")]));n.match2[1157] 
# n.match2[1157]<-n.match2[1157]+1;n.match2[1157]#adding the orig form itself

# #NOUN.ROOT[1160] #orig: 꽃잎을\구름빛/ncn+으로/jca", root 꽃잎을 'petal'
# n.match2[1160];n.match2[1160]<-c(sum(n.match2[1160], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="꽃잎")]));n.match2[1160] 

# NOUN.ROOT[1163] #orig 나를\혀/ncn" 
# n.match2[1163];n.match2[1163]<-c(sum(n.match2[1163], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="나")]));n.match2[1163]

# NOUN.ROOT[1169] #여름이었다/ncn", 여름 'summer'
# n.match2[1169];n.match2[1169]<-c(sum(n.match2[1169], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="여름")]));n.match2[1169]

# NOUN.ROOT[1182] #者바퀴벌레/ncn"
# n.match2[1182];n.match2[1182]<-c(kaist.n.freqs$V2[which(kaist.n.freqs$V1=="자")]);n.match2[1182] #need to change to 자 from 사람 ,note there are other meanings of 자  like 'ruler, foot (as in measurement)' also 'letter, character'
# #NOUN.ROOT[1183] #아이들혹은/ncn" , orig : "아이들,\\혹은,"
# n.match2[1183];n.match2[1183]<-c(sum(n.match2[1183], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="아이")]));n.match2[1183]

# NOUN.ROOT[1187] #orig: 세계에서\아무/ncn+도/jxc", root is 세계 'world'
# n.match2[1187];n.match2[1187]<-c(sum(n.match2[1187], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="세계")]));n.match2[1187]

# NOUN.ROOT[1193] #orig 물질\\거기 물질거/ncn+이/jp+기/etn", root 물질 'thing, matter, substance'
# n.match2[1193];n.match2[1193]<-c(sum(n.match2[1193], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="물질")]));n.match2[1193]

# NOUN.ROOT[1194] #그늘에서독하/ncn+이/jp+게/ecs, orig: 그늘에서\독하,  root "그늘" 'shade'
# n.match2[1194];n.match2[1194]<-c(sum(n.match2[1194], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="그늘")]));n.match2[1194]

# NOUN.ROOT[1199] #"반동이었다쓰러/ncn+이/jp+지/ef+고/jcr", orig: 반동이었다\쓰러, root 반동 'reaction'
# n.match2[1199];n.match2[1199]<-c(sum(n.match2[1199], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="반동")]));n.match2[1199]

# #NOUN.ROOT[1200]#것이냐털복숭/ncn+이/jcc" , orig = 것이냐\털복숭이
# n.match2[1200];n.match2[1200]<-c(sum(n.match2[1200], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="것")]));n.match2[1200]

# NOUN.ROOT[1204] #물방울소금처럼/ncn" orig: 물방울\\소금처럼, root 물방 'waterdrop'
# n.match2[1204];n.match2[1204]<-c(sum(n.match2[1204], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="물방")]));n.match2[1204]

# NOUN.ROOT[1205] #눈의한/ncn", orig "눈의\한" root 눈 'eyes'
# n.match2[1205];n.match2[1205]<-c(sum(n.match2[1205], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="눈")]));n.match2[1205]

# NOUN.ROOT[1206] #물은바람/ncn+에/jca, orig "물은\바람", root 'water'
# n.match2[1206];n.match2[1206]<-c(sum(n.match2[1206], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="물")]));n.match2[1206]

# NOUN.ROOT[1207] #우리졸/ncn", orig: 우리\졸  , root 우리 'us'
# n.match2[1207];n.match2[1207]<-c(sum(n.match2[1207], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="우리")]));n.match2[1207]

# NOUN.ROOT[1316] ("겨울나기") #겨울나/pvg+기/etn+에/jca"
# n.match2[1316];n.match2[1316]<-3;n.match2[1316] #manually found 2 examples of 겨울나기 in kaist.nouns3.txt + this example 

# NOUN.ROOT[1409]#<-c("눈썹") #눈썹，/ncn"  
# n.match2[1409];n.match2[1409]<-c(sum(n.match2[1409], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="눈썹")]));n.match2[1409]

# NOUN.ROOT[1436] #물，/ncn" 
# n.match2[1436];n.match2[1436]<-c(sum(n.match2[1436], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="물")]));n.match2[1436]

# NOUN.ROOT[1516] #목소/ncn" , orig 목소리가, root 목소리 'voice'
# n.match2[1516];n.match2[1516]<-c(sum(n.match2[1516], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="목소리")]));n.match2[1516]

# NOUN.ROOT[1520] #속삭이/pvg+ㅁ/etn"
# n.match2[1520];n.match2[1520]<-c(sum(n.match2[1520], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="속삭임")]));n.match2[1520]

# NOUN.ROOT[1578] #"아랫/ncn+날개/ncn+는/jxc" , "아랫날개" bottom wing   
# n.match2[1578];n.match2[1578]<-c(sum(n.match2[1578], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="아랫날개")]));n.match2[1578]

# NOUN.ROOT[1692] #길을/ncn"  
# n.match2[1692];n.match2[1692]<-c(sum(n.match2[1692], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="길")]));n.match2[1692]

# NOUN.ROOT[1780]#"펴/pvg+ㄴ/etm+지/nbn+로/jca" orig : 편지로  , root 편지 'letter'
# n.match2[1780];n.match2[1780]<-c(sum(n.match2[1780], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="편지")]));n.match2[1780]

# NOUN.ROOT[1788] #"아이로군……귀여운/ncn" 
# n.match2[1788];n.match2[1788]<-c(sum(n.match2[1788], kaist.n.freqs$V2[which(kaist.n.freqs$V1=="아이")]));n.match2[1788]
# NOUNFREQ<-n.match2

# ##fixing NOUN.ROOTS where there were compound nouns
# kaist.ncn.freqs<-read.table(file=("../perl/kaist.roots_ncn.txt"), quote="", sep="\t", comment.char="")
# kaist.ncn.freqs$V1<-as.vector(kaist.ncn.freqs$V1)
# ncn1<-exact.matches("^[^\\/]+", kaist.ncn.freqs$V1, pcre=TRUE, case.sens=FALSE)[[1]];head(ncn1);length(ncn1)
# ncn2<-exact.matches("ncn\\+[^\\/]+\\/ncn", kaist.ncn.freqs$V1, pcre=TRUE, case.sens=FALSE)[[1]];head(ncn2);length(ncn2)
# ncn2<-gsub("[\\d\\w\\s\\/\\+]", "", ncn2, perl=TRUE);head(ncn2);length(ncn2)
# ncnlemmas<-paste(ncn1, ncn2, sep="");head(ncnlemmas);length(ncnlemmas)
# kaist.ncn.freqs$V1<-ncnlemmas
# ncn.match<-c()
# for (j in 1:length(NOUN.ROOT2)){ 
		# if (sum(kaist.ncn.freqs$V1==NOUN.ROOT2[j])==0){
		# ncn.match <-c(ncn.match, 0)
		# } else {
			# ncn.match <-c(ncn.match, kaist.ncn.freqs$V2[which(kaist.ncn.freqs$V1==NOUN.ROOT2[j])])
			# }
	# }
# length(ncn.match); head(ncn.match)
# NOUN.ROOT2[which(ncn.match!=0)];ncn.match[which(ncn.match!=0)]
# NOUN.ROOT2[ncnposit[which(ncn.match==0)]]
# ncn_ncn[which(ncn.match==0)]
# NOUN_ANALYSIS[ncnposit[which(ncn.match==0)]]
# #will fix these two manually 국어/ncn+책/ncn	13 and the other is 음식/ncn+점/ncn 469
# ncn.match[477]<-13
# which(grepl("음식", NOUN.ROOT2))
# ncn.match[536]<-c(469+22)
# n.match3<-c(n.match2+ncn.match)
# which(n.match3==0)
# write.table(n.match3, file="../KAIST_n.match3.csv", quote=F, sep="\t", row.names=F, col.names=F)
n.match3 <-read.table(file="../KAIST_n.match3.csv", quote="", sep="\t", header=F, comment.char="")
n.match3<-n.match3$V1
NOUNFREQ<-n.match3

#continue with NSPEC FREQ 
#p(Adj|N) = p(Adj+N)/p(N) = f(Adj+N)/f(N)
#mod root noun root corpus tokens: kaist.modn.roots
#FIRST : GET p(Adj+N)
# MOD1NROOTS<-paste(MOD1.ROOT, NOUN.ROOT2, sep=" ")
# MOD2NROOTS<-paste(MOD2.ROOT, NOUN.ROOT2, sep=" ")
# mod1nmatches <-c()
# for (j in 1:length(MOD1NROOTS)){ 
		# mod1nmatches <-c(mod1nmatches, sum(grepl(MOD1NROOTS[j], kaist.modn.roots)))	
	# }
# head(mod1nmatches);length(mod1nmatches)#3911
# unique(MOD1NROOTS[which(mod1nmatches ==(max(mod1nmatches)))])
# hist(mod1nmatches)
# plot(sort(table(mod1nmatches)))
# mod2nmatches <-c()
# for (j in 1:length(MOD2NROOTS)){ 
		# mod2nmatches <-c(mod2nmatches, sum(grepl(MOD2NROOTS[j], kaist.modn.roots)))	
	# }
# head(mod2nmatches);length(mod2nmatches)#3911
# unique(MOD2NROOTS[which(mod2nmatches ==(max(mod2nmatches)))])
# hist(mod2nmatches)
# which(mod1nmatches==0);which(mod2nmatches==0)
# MOD2NROOTS[which(mod2nmatches==0)] #536 949
# MOD2NROOTS[949]#"차이 스웨터", orig 차이/ncn+나/ncn+식/xsncc (China-style)
# #fix this one manually
 # # kaist.files <- dir(".", full.names=TRUE) #
 # # chinastyle.matches <-c()
 # # chinastyle.matches2 <-c()
# # for (i in 1:length(kaist.files)){
  # # cat(i, "\n")
  # # curr.kaist.file<-scan(kaist.files[i], what=character(0), sep="\n", quiet=TRUE)
  # # # remove metadata lines
  # # curr.kaist.file<-gsub("<.*", "", curr.kaist.file, perl=TRUE)
  # # curr.chinastyle.matches <-exact.matches("차이나식.+", curr.kaist.file, pcre=TRUE, case.sens=FALSE)[[1]]
    # # curr.chinastyle.matches2 <-exact.matches("차이나식.+스웨터", curr.kaist.file, pcre=TRUE, case.sens=FALSE)[[1]]
 # # #this should get adjectives + noun roots
  # # # add to previous matches of adjective doublets
   # # chinastyle.matches <- c(chinastyle.matches, curr.chinastyle.matches)
   # # chinastyle.matches2 <- c(chinastyle.matches2, curr.chinastyle.matches2)
# # } 
# # head(chinastyle.matches);length(chinastyle.matches) #4
# # head(chinastyle.matches2);length(chinastyle.matches2)#1
# mod2nmatches[949]<-1
# #need to add counts for some mistagged things 
# mod1nmatches2<-mod1nmatches
# mod2nmatches2<-mod2nmatches
# #Fixing messed up NOUN.ROOTs 
# #n.match2[4]<-c("일") #일에대해 
# MOD1.ROOT[4];NOUN.ROOT2[4]
# mod1nmatches2[4];mod1nmatches2[4]<-c(sum(mod1nmatches2[4], length(which(kaist.modn.roots=="옳 일")))); mod1nmatches2[4]
# MOD2.ROOT[4];NOUN.ROOT2[4]
# mod2nmatches2[4];mod2nmatches2[4]<-c(sum(mod2nmatches2[4], length(which(kaist.modn.roots=="그르 일")))); mod2nmatches2[4]

# #NOUN.ROOT[83]<-c("저서")#저서우주찬미 
# MOD1.ROOT[83];NOUN.ROOT2[83]
# mod1nmatches2[83];mod1nmatches2[83]<-c(sum(mod1nmatches2[83], length(which(kaist.modn.roots=="작 저서")))); mod1nmatches2[83]
# MOD2.ROOT[83];NOUN.ROOT2[83]
# mod2nmatches2[83];mod2nmatches2[83]<-c(sum(mod2nmatches2[83], length(which(kaist.modn.roots=="아름답 저서")))); mod2nmatches2[83]

# #NOUN.ROOT[211]<-c("흐르")#흐르
# MOD1.ROOT[211];MOD2.ROOT[211]; NOUN.ROOT2[211]
# mod1nmatches2[211];mod1nmatches2[211]<-c(sum(mod1nmatches2[211], length(which(kaist.modn.roots=="파랗 흐름")))); mod1nmatches2[211]
# mod2nmatches2[211];mod2nmatches2[211]<-c(sum(mod2nmatches2[211], length(which(kaist.modn.roots=="많 흐름")))); mod2nmatches2[211]            

# #NOUN.ROOT[244]<-c("기쁨")#기쁘, 기쁘/paa+ㅁ/etn+을/jco
	# #need to add results for "기쁨" separately 
# MOD1.ROOT[244];MOD2.ROOT[244];NOUN.ROOT2[244]
# mod1nmatches2[244];mod1nmatches2[244]<-c(sum(mod1nmatches2[244], length(which(kaist.modn.roots=="크 기쁨")))); mod1nmatches2[244]
# mod2nmatches2[244];mod2nmatches2[244]<-c(sum(mod2nmatches2[244], length(which(kaist.modn.roots=="깊 기쁨")))); mod2nmatches2[244]  

# #NOUN.ROOT[315]<-c("생냉이")#생냉이가,
# MOD1.ROOT[315];MOD2.ROOT[315];NOUN.ROOT2[315]
# mod1nmatches2[315];mod1nmatches2[315]<-c(sum(mod1nmatches2[315], length(which(kaist.modn.roots=="검붉 생냉이")))); mod1nmatches2[315]
# mod2nmatches2[315];mod2nmatches2[315]<-c(sum(mod2nmatches2[315], length(which(kaist.modn.roots=="빳빳하 생냉이")))); mod2nmatches2[315]  

# #NOUN.ROOT[374]<-c("떨림")#떨리/pvg+ㅁ/etn+과/jct" 
# MOD1.ROOT[374];MOD2.ROOT[374];NOUN.ROOT2[374]
# mod1nmatches2[374];mod1nmatches2[374]<-c(sum(mod1nmatches2[374], length(which(kaist.modn.roots=="길 떨림")))); mod1nmatches2[374]
# mod2nmatches2[374];mod2nmatches2[374]<-c(sum(mod2nmatches2[374], length(which(kaist.modn.roots=="짧 떨림")))); mod2nmatches2[374]  

# #NOUN.ROOT[536] #"", "/sr+음식/ncn+점/ncn"
# #orig NOUN.ROOT is "" so MOD1NROOTS is "값싸 ", so has more matches than it should, need to replace it
# MOD1.ROOT[536];MOD2.ROOT[536];NOUN.ROOT2[536]
# mod1nmatches2[536];mod1nmatches2[536]<-c(length(which(kaist.modn.roots=="값싸 음식점"))); mod1nmatches2[536]
# mod2nmatches2[536];mod2nmatches2[536]<-c(length(which(kaist.modn.roots=="맛있 음식점"))); mod2nmatches2[536]  

# #NOUN.ROOT[900] #기쁘/paa+ㅁ/etn"
# MOD1.ROOT[900];MOD2.ROOT[900];NOUN.ROOT2[900]
# mod1nmatches2[900];mod1nmatches2[900]<-c(sum(mod1nmatches2[900], length(which(kaist.modn.roots=="뜨겁 기쁨")))); mod1nmatches2[900]
# mod2nmatches2[900];mod2nmatches2[900]<-c(sum(mod2nmatches2[900], length(which(kaist.modn.roots=="뭉클하 기쁨")))); mod2nmatches2[900]   
   
# #NOUN.ROOT[1037]#꽃들이오무라진/ncn"  add "꽃"
# MOD1.ROOT[1037];MOD2.ROOT[1037];NOUN.ROOT2[1037]
# mod1nmatches2[1037];mod1nmatches2[1037]<-c(sum(mod1nmatches2[1037], length(which(kaist.modn.roots=="빨갛 꽃")))); mod1nmatches2[1037]
# mod2nmatches2[1037];mod2nmatches2[1037]<-c(sum(mod2nmatches2[1037], length(which(kaist.modn.roots=="파랗 꽃")))); mod2nmatches2[1037]  

# #NOUN.ROOT[1090]#애를…/ncn+/sl"  
# MOD1.ROOT[1090];MOD2.ROOT[1090];NOUN.ROOT2[1090]
# mod1nmatches2[1090];mod1nmatches2[1090]<-c(sum(mod1nmatches2[1090], length(which(kaist.modn.roots=="이쁘 애")))); mod1nmatches2[1090]
# mod2nmatches2[1090];mod2nmatches2[1090]<-c(sum(mod2nmatches2[1090], length(which(kaist.modn.roots=="착하 애")))); mod2nmatches2[1090]   

# #NOUN.ROOT[1091] #아가씨였는데…/ncn+/sl"
# MOD1.ROOT[1091];MOD2.ROOT[1091];NOUN.ROOT2[1091]
# mod1nmatches2[1091];mod1nmatches2[1091]<-c(sum(mod1nmatches2[1091], length(which(kaist.modn.roots=="이쁘 아가씨")))); mod1nmatches2[1091]
# mod2nmatches2[1091];mod2nmatches2[1091]<-c(sum(mod2nmatches2[1091], length(which(kaist.modn.roots=="착하 아가씨")))); mod2nmatches2[1091]  

# NOUN.ROOT[1145] #비가\무너진,  noun is "비가", root "비" 'rain'
# MOD1.ROOT[1145];MOD2.ROOT[1145];NOUN.ROOT2[1145]
# mod1nmatches2[1145];mod1nmatches2[1145]<-c(sum(mod1nmatches2[1145], length(which(kaist.modn.roots=="희 비")))); mod1nmatches2[1145]
# mod2nmatches2[1145];mod2nmatches2[1145]<-c(sum(mod2nmatches2[1145], length(which(kaist.modn.roots=="엷 비")))); mod2nmatches2[1145]

# NOUN.ROOT[1149] #하늘에\새하얀 # 
# MOD1.ROOT[1149];MOD2.ROOT[1149];NOUN.ROOT2[1149]
# mod1nmatches2[1149];mod1nmatches2[1149]<-c(sum(mod1nmatches2[1149], length(which(kaist.modn.roots=="맑 하늘")))); mod1nmatches2[1149]
# mod2nmatches2[1149];mod2nmatches2[1149]<-c(sum(mod2nmatches2[1149], length(which(kaist.modn.roots=="푸르 하늘")))); mod2nmatches2[1149]  

# #NOUN.ROOT[1151] #새투명한, orig "새,\투명한" , root 새 'bird'
# MOD1.ROOT[1151];MOD2.ROOT[1151];NOUN.ROOT2[1151]
# mod1nmatches2[1151];mod1nmatches2[1151]<-c(sum(mod1nmatches2[1151], length(which(kaist.modn.roots=="귀엽 새")))); mod1nmatches2[1151]
# mod2nmatches2[1151];mod2nmatches2[1151]<-c(sum(mod2nmatches2[1151], length(which(kaist.modn.roots=="작 새")))); mod2nmatches2[1151]  

# #NOUN.ROOT[1153] #채찍이수없/ncn+이/jcc"  , orig = 채찍이\수없 , 채찍 'whip'
# MOD1.ROOT[1153];MOD2.ROOT[1153];NOUN.ROOT2[1153]
# mod1nmatches2[1153];mod1nmatches2[1153]<-c(sum(mod1nmatches2[1153], length(which(kaist.modn.roots=="매섭 채찍")))); mod1nmatches2[1153]
# mod2nmatches2[1153];mod2nmatches2[1153]<-c(sum(mod2nmatches2[1153], length(which(kaist.modn.roots=="길 채찍")))); mod2nmatches2[1153]  

# NOUN.ROOT[1155]<-c("하늘") #하늘을내보/ncn+이/jp+ㄴ다/ef", 하늘을내보 , orig = 하늘을\내보인다 , 하늘 'sky, heavens'
# MOD1.ROOT[1155];MOD2.ROOT[1155];NOUN.ROOT2[1155]
# mod1nmatches2[1155];mod1nmatches2[1155]<-c(sum(mod1nmatches2[1155], length(which(kaist.modn.roots=="작 하늘")))); mod1nmatches2[1155]
# mod2nmatches2[1155];mod2nmatches2[1155]<-c(sum(mod2nmatches2[1155], length(which(kaist.modn.roots=="깊 하늘")))); mod2nmatches2[1155] 

# NOUN.ROOT[1157] #시절\50원/ncn+을/jco" 
# MOD1.ROOT[1157];MOD2.ROOT[1157];NOUN.ROOT2[1157]
# mod1nmatches2[1157];mod1nmatches2[1157]<-c(sum(mod1nmatches2[1157], length(which(kaist.modn.roots=="벅차 시절")))); mod1nmatches2[1157]
# mod2nmatches2[1157];mod2nmatches2[1157]<-c(sum(mod2nmatches2[1157], length(which(kaist.modn.roots=="어리 시절")))); mod2nmatches2[1157]  

# #NOUN.ROOT[1160] #orig: 꽃잎을\구름빛/ncn+으로/jca", root 꽃잎을 'petal'
# MOD1.ROOT[1160];MOD2.ROOT[1160];NOUN.ROOT2[1160]
# mod1nmatches2[1160];mod1nmatches2[1160]<-c(sum(mod1nmatches2[1160], length(which(kaist.modn.roots=="작 꽃잎")))); mod1nmatches2[1160]
# mod2nmatches2[1160];mod2nmatches2[1160]<-c(sum(mod2nmatches2[1160], length(which(kaist.modn.roots=="여리 꽃잎")))); mod2nmatches2[1160]  

# NOUN.ROOT[1163] #orig 나를\혀/ncn" 
# MOD1.ROOT[1163];MOD2.ROOT[1163];NOUN.ROOT2[1163]
# mod1nmatches2[1163];mod1nmatches2[1163]<-c(sum(mod1nmatches2[1163], length(which(kaist.modn.roots=="하찮 나")))); mod1nmatches2[1163]
# mod2nmatches2[1163];mod2nmatches2[1163]<-c(sum(mod2nmatches2[1163], length(which(kaist.modn.roots=="물렁하 나")))); mod2nmatches2[1163]  

# NOUN.ROOT[1169] #여름이었다/ncn", 여름 'summer'
# MOD1.ROOT[1169];MOD2.ROOT[1169];NOUN.ROOT2[1169]
# mod1nmatches2[1169];mod1nmatches2[1169]<-c(sum(mod1nmatches2[1169], length(which(kaist.modn.roots=="무덥 여름")))); mod1nmatches2[1169]
# mod2nmatches2[1169];mod2nmatches2[1169]<-c(sum(mod2nmatches2[1169], length(which(kaist.modn.roots=="길 여름")))); mod2nmatches2[1169]  

# NOUN.ROOT[1182] #者바퀴벌레/ncn"
# MOD1.ROOT[1182];MOD2.ROOT[1182];NOUN.ROOT2[1182]
# mod1nmatches2[1182];mod1nmatches2[1182]<-c(sum(mod1nmatches2[1182], length(which(kaist.modn.roots=="배고프 자")))); mod1nmatches2[1182]
# mod2nmatches2[1182];mod2nmatches2[1182]<-c(sum(mod2nmatches2[1182], length(which(kaist.modn.roots=="외롭 자")))); mod2nmatches2[1182]  

# #NOUN.ROOT[1183] #아이들혹은/ncn" , orig : "아이들,\\혹은,"
# MOD1.ROOT[1183];MOD2.ROOT[1183];NOUN.ROOT2[1183]
# mod1nmatches2[1183];mod1nmatches2[1183]<-c(sum(mod1nmatches2[1183], length(which(kaist.modn.roots=="많 아이")))); mod1nmatches2[1183]
# mod2nmatches2[1183];mod2nmatches2[1183]<-c(sum(mod2nmatches2[1183], length(which(kaist.modn.roots=="어리 아이")))); mod2nmatches2[1183]  

# NOUN.ROOT[1187] #orig: 세계에서\아무/ncn+도/jxc", root is 세계 'world'
# MOD1.ROOT[1187];MOD2.ROOT[1187];NOUN.ROOT2[1187]
# mod1nmatches2[1187];mod1nmatches2[1187]<-c(sum(mod1nmatches2[1187], length(which(kaist.modn.roots=="어둡 세계")))); mod1nmatches2[1187]
# mod2nmatches2[1187];mod2nmatches2[1187]<-c(sum(mod2nmatches2[1187], length(which(kaist.modn.roots=="축축하 세계")))); mod2nmatches2[1187]  

# NOUN.ROOT[1193] #orig 물질\\거기 물질거/ncn+이/jp+기/etn", root 물질 'thing, matter, substance'
# MOD1.ROOT[1193];MOD2.ROOT[1193];NOUN.ROOT2[1193]
# mod1nmatches2[1193];mod1nmatches2[1193]<-c(sum(mod1nmatches2[1193], length(which(kaist.modn.roots=="미끄럽 물질")))); mod1nmatches2[1193]
# mod2nmatches2[1193];mod2nmatches2[1193]<-c(sum(mod2nmatches2[1193], length(which(kaist.modn.roots=="차갑 물질")))); mod2nmatches2[1193]  

# NOUN.ROOT[1194] #그늘에서독하/ncn+이/jp+게/ecs, orig: 그늘에서\독하,  root "그늘" 'shade'
# MOD1.ROOT[1194];MOD2.ROOT[1194];NOUN.ROOT2[1194]
# mod1nmatches2[1194];mod1nmatches2[1194]<-c(sum(mod1nmatches2[1194], length(which(kaist.modn.roots=="비좁 그늘")))); mod1nmatches2[1194]
# mod2nmatches2[1194];mod2nmatches2[1194]<-c(sum(mod2nmatches2[1194], length(which(kaist.modn.roots=="컴컴하 그늘")))); mod2nmatches2[1194]  

# NOUN.ROOT[1199] #"반동이었다쓰러/ncn+이/jp+지/ef+고/jcr", orig: 반동이었다\쓰러, root 반동 'reaction'
# MOD1.ROOT[1199];MOD2.ROOT[1199];NOUN.ROOT2[1199]
# mod1nmatches2[1199];mod1nmatches2[1199]<-c(sum(mod1nmatches2[1199], length(which(kaist.modn.roots=="아름답 반동")))); mod1nmatches2[1199]
# mod2nmatches2[1199];mod2nmatches2[1199]<-c(sum(mod2nmatches2[1199], length(which(kaist.modn.roots=="눈물겹 반동")))); mod2nmatches2[1199] 

# #NOUN.ROOT[1200]#것이냐털복숭/ncn+이/jcc" , orig = 것이냐\털복숭이
# MOD1.ROOT[1200];MOD2.ROOT[1200];NOUN.ROOT2[1200]
# mod1nmatches2[1200];mod1nmatches2[1200]<-c(sum(mod1nmatches2[1200], length(which(kaist.modn.roots=="부끄럽 것")))); mod1nmatches2[1200]
# mod2nmatches2[1200];mod2nmatches2[1200]<-c(sum(mod2nmatches2[1200], length(which(kaist.modn.roots=="덧없 것")))); mod2nmatches2[1200]  

# NOUN.ROOT[1204] #물방울소금처럼/ncn" orig: 물방울\\소금처럼, root 물방 'waterdrop'
# MOD1.ROOT[1204];MOD2.ROOT[1204];NOUN.ROOT2[1204]
# mod1nmatches2[1204];mod1nmatches2[1204]<-c(sum(mod1nmatches2[1204], length(which(kaist.modn.roots=="작 물방")))); mod1nmatches2[1204]
# mod2nmatches2[1204];mod2nmatches2[1204]<-c(sum(mod2nmatches2[1204], length(which(kaist.modn.roots=="푸르 물방")))); mod2nmatches2[1204]  

# NOUN.ROOT[1205] #눈의한/ncn", orig "눈의\한" root 눈 'eyes'
# MOD1.ROOT[1205];MOD2.ROOT[1205];NOUN.ROOT2[1205]
# mod1nmatches2[1205];mod1nmatches2[1205]<-c(sum(mod1nmatches2[1205], length(which(kaist.modn.roots=="깊 눈")))); mod1nmatches2[1205]
# mod2nmatches2[1205];mod2nmatches2[1205]<-c(sum(mod2nmatches2[1205], length(which(kaist.modn.roots=="푸르 눈")))); mod2nmatches2[1205]  

# NOUN.ROOT[1206] #물은바람/ncn+에/jca, orig "물은\바람", root 'water'
# MOD1.ROOT[1206];MOD2.ROOT[1206];NOUN.ROOT2[1206]
# mod1nmatches2[1206];mod1nmatches2[1206]<-c(sum(mod1nmatches2[1206], length(which(kaist.modn.roots=="높 물")))); mod1nmatches2[1206]
# mod2nmatches2[1206];mod2nmatches2[1206]<-c(sum(mod2nmatches2[1206], length(which(kaist.modn.roots=="너르 물")))); mod2nmatches2[1206]  

# NOUN.ROOT[1207] #우리졸/ncn", orig: 우리\졸  , root 우리 'us'
# MOD1.ROOT[1207];MOD2.ROOT[1207];NOUN.ROOT2[1207]
# mod1nmatches2[1207];mod1nmatches2[1207]<-c(sum(mod1nmatches2[1207], length(which(kaist.modn.roots=="잘나 우리")))); mod1nmatches2[1207]
# mod2nmatches2[1207];mod2nmatches2[1207]<-c(sum(mod2nmatches2[1207], length(which(kaist.modn.roots=="뻣뻣하 우리")))); mod2nmatches2[1207]  

# NOUN.ROOT[1316] #("겨울나기") #겨울나/pvg+기/etn+에/jca"
# MOD1.ROOT[1316];MOD2.ROOT[1316];NOUN.ROOT2[1316]
# mod1nmatches2[1316];mod1nmatches2[1316]<-c(sum(mod1nmatches2[1316], length(which(kaist.modn.roots=="춥 겨울나기")))); mod1nmatches2[1316]
# mod2nmatches2[1316];mod2nmatches2[1316]<-c(sum(mod2nmatches2[1316], length(which(kaist.modn.roots=="배고프 겨울나기")))); mod2nmatches2[1316]  

# NOUN.ROOT[1409]#<-c("눈썹") #눈썹，/ncn"  
# MOD1.ROOT[1409];MOD2.ROOT[1409];NOUN.ROOT2[1409]
# mod1nmatches2[1409];mod1nmatches2[1409]<-c(sum(mod1nmatches2[1409], length(which(kaist.modn.roots=="가늘 눈썹")))); mod1nmatches2[1409]
# mod2nmatches2[1409];mod2nmatches2[1409]<-c(sum(mod2nmatches2[1409], length(which(kaist.modn.roots=="길 눈썹")))); mod2nmatches2[1409]  

# NOUN.ROOT[1436] #물，/ncn" 
# MOD1.ROOT[1436];MOD2.ROOT[1436];NOUN.ROOT2[1436]
# mod1nmatches2[1436];mod1nmatches2[1436]<-c(sum(mod1nmatches2[1436], length(which(kaist.modn.roots=="맑 물")))); mod1nmatches2[1436]
# mod2nmatches2[1436];mod2nmatches2[1436]<-c(sum(mod2nmatches2[1436], length(which(kaist.modn.roots=="차갑 물")))); mod2nmatches2[1436]  

# NOUN.ROOT[1516] #목소/ncn" , orig 목소리가, root 목소리 'voice'
# MOD1.ROOT[1516];MOD2.ROOT[1516];NOUN.ROOT2[1516]
# mod1nmatches2[1516];mod1nmatches2[1516]<-c(sum(mod1nmatches2[1516], length(which(kaist.modn.roots=="높 목소리")))); mod1nmatches2[1516]
# mod2nmatches2[1516];mod2nmatches2[1516]<-c(sum(mod2nmatches2[1516], length(which(kaist.modn.roots=="날카롭 목소리")))); mod2nmatches2[1516]

# NOUN.ROOT[1520] #속삭이/pvg+ㅁ/etn"
# MOD1.ROOT[1520];MOD2.ROOT[1520];NOUN.ROOT2[1520]
# mod1nmatches2[1520];mod1nmatches2[1520]<-c(sum(mod1nmatches2[1520], length(which(kaist.modn.roots=="부드럽 속삭임")))); mod1nmatches2[1520]
# mod2nmatches2[1520];mod2nmatches2[1520]<-c(sum(mod2nmatches2[1520], length(which(kaist.modn.roots=="나긋나긋하 속삭임")))); mod2nmatches2[1520]  
 
# NOUN.ROOT[1692] #길을/ncn"  
# MOD1.ROOT[1692];MOD2.ROOT[1692];NOUN.ROOT2[1692]
# mod1nmatches2[1692];mod1nmatches2[1692]<-c(sum(mod1nmatches2[1692], length(which(kaist.modn.roots=="넓 길")))); mod1nmatches2[1692]
# mod2nmatches2[1692];mod2nmatches2[1692]<-c(sum(mod2nmatches2[1692], length(which(kaist.modn.roots=="딱딱하 길")))); mod2nmatches2[1692]  

# NOUN.ROOT[1780]#"펴/pvg+ㄴ/etm+지/nbn+로/jca" orig : 편지로  , root 편지 'letter'
# MOD1.ROOT[1780];MOD2.ROOT[1780];NOUN.ROOT2[1780]
# mod1nmatches2[1780];mod1nmatches2[1780]<-c(sum(mod1nmatches2[1780], length(which(kaist.modn.roots=="길 편지")))); mod1nmatches2[1780]
# mod2nmatches2[1780];mod2nmatches2[1780]<-c(sum(mod2nmatches2[1780], length(which(kaist.modn.roots=="펴 편지")))); mod2nmatches2[1780]  

# NOUN.ROOT[1788] #"아이로군……귀여운/ncn" 
# MOD1.ROOT[1788];MOD2.ROOT[1788];NOUN.ROOT2[1788]
# mod1nmatches2[1788];mod1nmatches2[1788]<-c(sum(mod1nmatches2[1788], length(which(kaist.modn.roots=="귀엽 아이")))); mod1nmatches2[1788]
# mod2nmatches2[1788];mod2nmatches2[1788]<-c(sum(mod2nmatches2[1788], length(which(kaist.modn.roots=="잘생기 아이")))); mod2nmatches2[1788]  

# #Fix the rest
# MOD1.ROOT[which(mod2nmatches==0)];MOD2.ROOT[which(mod2nmatches==0)];NOUN.ROOT2[which(mod2nmatches==0)]
# which(mod2nmatches2==0) #none
# MOD2_ANALYSIS[which(mod2nmatches2==0)] #차이/ncn+나/ncn+식/xsncc
# MOD2NROOTS[which(mod2nmatches2==0)]
# write.table(mod1nmatches2, file="../KAIST_mod1nmatches2.csv", quote=F, sep="\t", row.names=F, col.names=F)
mod1nmatches2<-read.table(file="../KAIST_mod1nmatches2.csv", quote="", sep="\t", header=F)
mod1nmatches2<-mod1nmatches2$V1
#write.table(mod2nmatches2, file="../KAIST_mod2nmatches2.csv", quote=F, sep="\t", row.names=F, col.names=F)
mod2nmatches2<-read.table(file="../KAIST_mod2nmatches2.csv", quote="", sep="\t", header=F)
mod2nmatches2<-mod2nmatches2$V1
NSPEC1<-c(mod1nmatches2/n.match3) #one word, 아프, has value higher than 1
NSPEC2<-c(mod2nmatches2/n.match3)
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
# sum(grepl("^\\W+\\/n", MOD1_ANALYSIS, perl=TRUE)) #0
# sum(grepl("^\\W+\\/p", MOD1_ANALYSIS, perl=TRUE)) #711
# sum(grepl("^\\W+\\/n", MOD2_ANALYSIS, perl=TRUE)) #8
# sum(grepl("^\\W+\\/p", MOD2_ANALYSIS, perl=TRUE)) #703
# NOMCHAR1<-grepl("^\\W+/n", MOD1_ANALYSIS, perl=TRUE)
# NOMCHAR2<-grepl("^\\W+/n", MOD2_ANALYSIS, perl=TRUE)
# NOMCHAR1[which(NOMCHAR1=="FALSE")]<-c("adj")
# NOMCHAR1[which(NOMCHAR1=="TRUE")]<-c("n")
# NOMCHAR2[which(NOMCHAR2=="FALSE")]<-c("adj")
# NOMCHAR2[which(NOMCHAR2=="TRUE")]<-c("n")
# table(NOMCHAR1);table(NOMCHAR2)
# #updated to Wona_Coding_final.csv, there is only one n match in the whole corpus, so should investigate this further but this measure is useless for now
# length(grep("\\/n", kaist.analyzed$MOD2_ANALYSIS, perl=TRUE)) #202
# length(grep("\\/p", kaist.analyzed$MOD2_ANALYSIS, perl=TRUE)) #3791
# #NOMCHAR (END)
# SEMOPEN1<-semopen$SEMOPEN1
# SEMOPEN2<-semopen$SEMOPEN2
# INDCOMP1<-indcomp$INDCOMP1
# INDCOMP2<-indcomp$INDCOMP2
# FILE<-kaist.mod.mod.n$FILE
#AFFLOAD1<-MOD1_AFFLOAD
#AFFLOAD2<-MOD2_AFFLOAD

#SUBJOBJ (START) 
#Ones to fix
# MOD1_SUBJOBJ[72] <-6#CASE=149, 곧고 긴 움길을,  곧고, 'straight' MOD1_SUBJOBJ = 5, should be 6
# MOD2_SUBJOBJ[162]<-5#"high-blue" should be 5, was accidentally put as 12-5
# MOD2_SUBJOBJ[224]<-7#CASE = 433, 낡고 오래된 트럭들이었으며 'rugged and old truck', old is 8, should be 7
# MOD1_SUBJOBJ[248]<-11#CASE=480, 밝고 더운 기운이 'bright and hot atmosphere', bright is 12, should be 11
# MOD1_SUBJOBJ[493]<-7 #CASE= 930, 젊고 예쁜 여성이라면 'young and pretty lady', young is 12, should be 7
# MOD2_SUBJOBJ[584]<-5#CASE=1155, 빨갛고 노란 인조,  노란 'yellow' MOD2_SUBJOBJ = 2, should be 5
# MOD1_SUBJOBJ[666]<-12#CASE=1336 무섭고 지겨운 개강을 'scary and boring new semester', scary and boring are both 8, should be 12
# MOD2_SUBJOBJ[666]<-12
# #not added: #CASE=1607 춥고 배고픈 산적들이 'cold and hungry bandit', hungry SUBJOBJ was put into DEF accidentally (watch this one, might not have been added in the first place), should be 8
# MOD1_SUBJOBJ[209]<-12#CASE=412 싫고 지겨운 일이라고 'to be disliked and boring job', both are 8 , should be 12
# MOD2_SUBJOBJ[209]<-12
# MOD1_SUBJOBJ[1567]<-12#CASE=3338 귀찮고 성가신 일인 'annoying and bothersome work', both are 8, should be 12
# MOD2_SUBJOBJ[1567]<-12
# # not added MOD1_SUBJOBJ[72]#CASE=2345 즐거운 서글픈 표정을  happy and sad facial-expression, AFFLOAD was put into SUBJOBJ accidentally (watch this one, might not have been added in the first place), SUBJOBJ should be 8
# MOD2_SUBJOBJ[1253]<-12#CASE=2489, 새로운 짜릿한 만남 'new and exciting meeting', 'exciting' is 8, should be 12
# MOD1_SUBJOBJ[1397]<-11#CASE=2870 #짧거나 긴 것이 ,짧거나 'short' was accidentally put as 1 instead of 11 (MOD1_SUBJOBJ)
# MOD1_SUBJOBJ[1545]<-7#CASE=3286, 오래된 낡은 차를 'old and rugged car', 'old' is 11, should be 7
# MOD1_SUBJOBJ[1791]<-10#CASE=3849, 질기고 맛없는 고기밖에, 'tough and not-tasty meat' , tough is 11, should be 10 (matches other example of this)
# #a few are blank that probably shouldn't be
# MOD1_SUBJOBJ<-factor(MOD1_SUBJOBJ)
# MOD2_SUBJOBJ<-factor(MOD2_SUBJOBJ)
# SUBJOBJ1<-as.vector(MOD1_SUBJOBJ)
# SUBJOBJ2<-as.vector(MOD2_SUBJOBJ)
# SUBJOBJ1<-as.numeric(SUBJOBJ1)
# SUBJOBJ2<-as.numeric(SUBJOBJ2)
#changing some dark color ones from 5 to 11, in Wulff's discussion of Dixon's category ordering she uses "light blue paper" as an example of preferred ordering for category 5 physical property, so that would indicate she would agree in not considering "light" in "light blue paper" to be a color but rather a physical property
# darkpos<-(grep("dark", MOD1_DEF, perl=T))
# PHRASE_DEF[darkpos]
# SUBJOBJ1[darkpos]
# PHRASE_DEF[darkpos[SUBJOBJ1[darkpos]]]
# SUBJOBJ1[darkpos[SUBJOBJ1[darkpos]==5]]<-11
# darkpos<-(grep("dark", MOD2_DEF, perl=T))
# PHRASE_DEF[darkpos];SUBJOBJ2[darkpos]
# PHRASE_DEF[darkpos[SUBJOBJ2[darkpos]==5]]
# SUBJOBJ2[darkpos[SUBJOBJ2[darkpos]==5]]<-11
# #"light" is always 11 or 12
# brightpos<-(grep("bright", MOD1_DEF, perl=T))
# PHRASE_DEF[brightpos];SUBJOBJ1[brightpos]
# PHRASE_DEF[brightpos[SUBJOBJ1[brightpos]==5]] # only 1 example
# SUBJOBJ1[brightpos[SUBJOBJ1[brightpos]==5]]<-11
# #changing new and old to Age, in Wulff's discussion of Dixon's category ordering she indicates "new" as in "important new publication" as Age category
# newpos<-(grep("new", MOD1_DEF, perl=T))
# PHRASE_DEF[newpos]; SUBJOBJ1[newpos]
# PHRASE_DEF[newpos[SUBJOBJ1[newpos]!=7]] # two were coded with 7 
# SUBJOBJ1[newpos[SUBJOBJ1[newpos]!=7]]<-7 
# newpos <-(grep("new", MOD2_DEF, perl=T))
# PHRASE_DEF[newpos]; SUBJOBJ2[newpos]
# PHRASE_DEF[newpos[SUBJOBJ2[newpos]==7]] # only 1 example
# SUBJOBJ2[newpos[SUBJOBJ2[newpos]!=7]]<-7
# #23/24 of "new" in SUBJOBJ2 position are from "unique and new fact" , could think about redoing the model with phrases with a certain threshold of repetitions reduced to one token
# oldpos<-(grep("^old", MOD1_DEF, perl=T))
# PHRASE_DEF[oldpos]; SUBJOBJ1[oldpos]
# PHRASE_DEF[oldpos[SUBJOBJ1[oldpos]!=7]] 
# SUBJOBJ1[oldpos[SUBJOBJ1[oldpos]!=7]]<-7
# oldpos <-which(MOD2_DEF=="old")
# PHRASE_DEF[oldpos]; SUBJOBJ2[oldpos]
# PHRASE_DEF[oldpos[SUBJOBJ2[oldpos]!=7]] 
# SUBJOBJ2[oldpos[SUBJOBJ2[oldpos]!=7]]<-7
# #MOD1 "young" are all 7 
# youngpos <-which(MOD2_DEF=="young")
# PHRASE_DEF[youngpos]; SUBJOBJ2[youngpos]
# PHRASE_DEF[youngpos[SUBJOBJ2[youngpos]!=7]] #뜨거운 젊은 피를 hot and young blood 
# SUBJOBJ2[youngpos[SUBJOBJ2[youngpos]!=7]]<-7
# #SUBJOBJ (END)

#kaistresults<-data.frame(FILE, MOD1, RT.LENGTH1, SEMOPEN1, INDCOMP1, SUBJOBJ1, MOD1_AFFLOAD, NSPEC1, FREQ1, MOD2, RT.LENGTH2, SEMOPEN2, INDCOMP2, SUBJOBJ2, MOD2_AFFLOAD, NSPEC2, FREQ2, NOUN, NOUNFREQ)
#write.table(kaistresults, file="../KAIST_finalresults.csv", quote=F, sep="\t", row.names=F, col.names=T)
#kaistresults<-read.table(file="../KAIST_finalresults.csv", quote="", sep="\t", header=T)
#attach(kaistresults)
# PHRASE<-paste(MOD1, MOD2, NOUN, sep=" ")
#PHRASE_DEF<-Definition
#kaist.everything<-data.frame(FILE, PHRASE, PHRASE_DEF, MOD1, MOD1.ROOT, MOD1_ANALYSIS, MOD1_DEF, CONNECT, RT.LENGTH1, SEMOPEN1, INDCOMP1, SUBJOBJ1, MOD1_AFFLOAD, NSPEC1, FREQ1, MOD2, MOD2.ROOT, MOD2_ANALYSIS, MOD2_DEF, RT.LENGTH2, SEMOPEN2, INDCOMP2, SUBJOBJ2, MOD2_AFFLOAD, NSPEC2, FREQ2, NOUN, NOUN.ROOT, NOUN_ANALYSIS, NOUNFREQ)
#write.table(kaist.everything, file="../KAIST_everything.csv", quote=F, sep="\t", row.names=F, col.names=T)
kaisteverything<-read.table(file="../KAIST_everything.csv", quote="", sep="\t", header=T)
attach(kaisteverything)

length(unique(FILE)) #780
x<-hist(table(FILE)) #719 files appeared between 0 and 5 times
719/780 #92.2%
file.freq<-c()
for (j in 1:length(FILE)){ 
		file.freq <-c(file.freq, length(which(FILE==FILE[j])))	
	}
 head(file.freq);length(file.freq)#1810
 unique(FILE[which(file.freq ==(max(file.freq)))])
uniq.files<-data.frame(FILE, file.freq)
uniq.files<-unique(uniq.files)
sum(uniq.files$file.freq[which(uniq.files$file.freq<=5)]) # = 1241 = number of data points represented by files with freq between 1 and 5 
1241/1810 #68.6%

#Let's see what the type/token ratio looks like for phrases 
ROOTPHRASE<-paste(MOD1.ROOT, MOD2.ROOT, NOUN.ROOT, sep=" ")
length(ROOTPHRASE);length(unique(ROOTPHRASE)) #1625 types for 1810 tokens
tail(sort(table(ROOTPHRASE)))
#새롭 낯설 경험  18 "new and unfamiliar 
#옳 그르 것  15
#좋 나쁘 것  12
#작 앳되 얼굴 4
#새롭 낯설 사람들  4
#비싸 좋 것  4
hist(table(ROOTPHRASE))

##REGRESSION ANALYSIS
#LENGTH predicted: longer modifier will be closer to noun
#LENGTH 1 < LENGTH 2
RTLENGTHDIFF<-c(RT.LENGTH2-RT.LENGTH1, RT.LENGTH1-RT.LENGTH2);length(RTLENGTHDIFF)
#No NOMCHAR
#SEMOPEN predicted: less semantically open (modifying less noun types) closer to n
#SEMOPEN 1 > SEMOPEN 2 NOUN
SEMOPENDIFF<-c(SEMOPEN1-SEMOPEN2, SEMOPEN2-SEMOPEN1)
#INDCOMP predicted: more independent from comparison (low value of INDCOMP), closer to n
#INDCOMP 1 > INDCOMP 2 NOUN
INDCOMPDIFF<-c(INDCOMP1-INDCOMP2, INDCOMP2-INDCOMP1)
#SUBOBJ predicted: more objective (lower value) will be closer to noun
#SUBJOBJ1 > SUBJOBJ2
SUBJOBJDIFF<-c(SUBJOBJ1-SUBJOBJ2, SUBJOBJ2-SUBJOBJ1)
#AFFLOAD1 > AFFLOAD2
AFFLOADDIFF<-c(MOD1_AFFLOAD-MOD2_AFFLOAD, MOD2_AFFLOAD-MOD1_AFFLOAD)
#NSPECFREQ predicted: high value is high association, highly associated adj is predicted to appear closer to the noun
#p(Adj|N) = p(Adj+N)/p(N) = f(Adj+N)/f(N) : conditional probability of knowing the adjectives given knowledge of the noun
#NSPEC1 < NSPEC 2
NSPECDIFF<-c(NSPEC2-NSPEC1,NSPEC1-NSPEC2)
#GENFREQ predicted: more frequent adj precede less frequent ones , so lower value would be closer to noun
#FREQ1 > FREQ2  (logged values of GENFREQ)
GENFREQDIFF<-c(FREQ1-FREQ2, FREQ2-FREQ1)
dep<-c(rep("predict", 1810), rep("nonpredict", 1810))
dep<-as.factor(dep)
#Variable for connecting grammar
# CONNECT<-c(rep("은/는", 1810))
# CONNECT[grep("고$", MOD1, perl=T)]<-c("고")
# CONNECT[grep("거나$", MOD1, perl=T)]<-c("거나")
table(CONNECT)
# 거나    고 은/는 
#  44  1465   301 
CONNECT[mismatch]
table(CONNECT[mismatch])
# 거나    고 은/는 
#    6   111    20 
#The proportion of 거나 that is in mismatch  is 0.1363636 , about twice that of the other two (0.07576792 for 고  and 0.06644518 for 은/는 ), so there could be an effect there? It's still pretty low tho
# #CONNECT<-(CONNECT, CONNECT) #?? the problem is what should be the value for non-preferred, also the
CONNECT2 <-CONNECT
CONNECT2 <- factor(CONNECT2, levels = c(levels(CONNECT2), "etm", "ecc"))
CONNECT2[grep("은", CONNECT)]<-c("etm")
CONNECT2[grep("거나", CONNECT)]<-c("ecc")
CONNECT2[grep("고", CONNECT)]<-c("ecc")
CONNECT2 <- factor(CONNECT2, levels = c("etm", "ecc"))
table(CONNECT2)
CONNECTDIFF<-append(CONNECT2, CONNECT2)
CONNECTDIFF <-factor(CONNECTDIFF, levels=1:nlevels(CONNECT2), labels=levels(CONNECT2))
table(CONNECTDIFF)


#Make version with shorter names for plot
RLdiff<-RTLENGTHDIFF
RLdiff<-as.vector(RLdiff)
SemOpdiff<-SEMOPENDIFF
SemOpdiff<-as.vector(SemOpdiff)
ICdiff<-INDCOMPDIFF
ICdiff<-as.vector(ICdiff)
SubObjdiff<-SUBJOBJDIFF
SubObjdiff<-as.vector(SubObjdiff)
ALdiff<-AFFLOADDIFF
ALdiff<-as.vector(ALdiff)
NSpFdiff<-NSPECDIFF
NSpFdiff<-as.vector(NSpFdiff)
GFdiff<-GENFREQDIFF
GFdiff<-as.vector(GFdiff)
CON<-CONNECTDIFF
ORDER<-dep
model2.glm6<-glm(ORDER~(RLdiff+ SemOpdiff + ICdiff + SubObjdiff + ALdiff + NSpFdiff + GFdiff + CON + RLdiff: CON + ICdiff: CON + NSpFdiff:CON + SemOpdiff: CON + SubObjdiff: CON), family=binomial)
summary(model2.glm6)
# Coefficients:
                    # Estimate Std. Error z value Pr(>|z|)    
# (Intercept)        1.814e-15  1.266e-01   0.000  1.00000    
# RLdiff            -7.001e-01  1.341e-01  -5.221 1.78e-07 ***
# SemOpdiff          5.359e+00  7.983e-01   6.713 1.91e-11 ***
# ICdiff             2.662e+01  3.173e+00   8.391  < 2e-16 ***
# SubObjdiff         3.266e-02  4.066e-02   0.803  0.42176    
# ALdiff             3.103e-01  6.963e-02   4.456 8.34e-06 ***
# NSpFdiff           9.444e+01  1.268e+01   7.447 9.52e-14 ***
# GFdiff             2.166e-01  3.159e-02   6.856 7.06e-12 ***
# CONecc            -1.170e-15  1.331e-01   0.000  1.00000    
# RLdiff:CONecc      1.366e+00  1.424e-01   9.592  < 2e-16 ***
# ICdiff:CONecc     -2.432e+01  3.341e+00  -7.278 3.40e-13 ***
# NSpFdiff:CONecc   -7.500e+01  1.280e+01  -5.860 4.63e-09 ***
# SemOpdiff:CONecc  -3.691e+00  8.012e-01  -4.607 4.09e-06 ***
# SubObjdiff:CONecc -1.366e-01  4.373e-02  -3.123  0.00179 ** 

library(effects)
par(cex.main=0.1, cex.axis=0.5)
plot(allEffects(model2.glm6), ask=FALSE, grid=TRUE, ylim=c(0,1), rescale.axis=FALSE, cex=1.0) 
#plotting effects separately
plot(effect("NSpFdiff:CON", model2.glm6))
plot(effect("SubObjdiff:CON", model2.glm6))
#looking at max of each variable: 
max(RTLENGTHDIFF);max(SEMOPENDIFF);max(INDCOMPDIFF);max(SUBJOBJDIFF);max(AFFLOADDIFF);max(NSPECDIFF);max(GENFREQDIFF)
#Classification accuracy
somers2(binomial()$linkinv(fitted(model2.glm6)), as.numeric(dep)-1) 
#C value = 0.7949879
# to get Nagelkerke's R2 and C
library(rms)
lrm(formula(model2.glm6)) #R2 = 0.367
actual<-data.frame(RLdiff, SemOpdiff, ICdiff, SubObjdiff, ALdiff, NSpFdiff, GFdiff, CONdiff)
classifications<-predict(model2.glm6, type="response")
classif2<-ifelse(classifications>=0.5, "predict", "nonpredict")
evaluation<-table(classif2, dep)
addmargins(evaluation);sum(diag(evaluation))/sum(evaluation) #70% correct

#which ones are not predicted correctly? 
mismatch<-c(unique(which(classif2!=dep)))
mismatch<-c(mismatch[1:540])
MOD1[mismatch]
length(which(RTLENGTHDIFF[mismatch]<0)) #240
SUBJOBJ1[mismatch]; SUBJOBJ2[mismatch]
SUBJOBJDIFF[mismatch]
length(which(SUBJOBJDIFF[mismatch]<0)) #111
length(sort(union((which(SUBJOBJDIFF[mismatch]<0)), (which(RTLENGTHDIFF[mismatch]<0)))))
#so 279 of the 540 mismatches involve one or both of these 
length(which(SUBJOBJDIFF[1:1810]<0)) #432 

 #SEMOPEN thoughts
 length(which(kaist.modn.modroot=="크"))
length(grep("크\\s", kaist.modn.unique))
length(which(kaist.modn.modroot.unique=="크"))
 #it seems weird to penalize "크" (which ends up having a really low SEMOPEN value) for having many repeated items with the same noun, that is the only thing that relativizing SEMOPEN by the freq of the adjective does is penalize ones that have many of the same noun repeated, however doesn't the fact of them appearing many diff noun types mean that they are semantically open in itself? it shouldn't matter if some of those noun types end up being repeated a lot, the point is that if the adj appears with many types in the corpus, it is able to be used with many diff noun types, it seems to weaken the measure to include information about how many of those are common repeated phrases
 #in fact: 
semopenboth<-c(SEMOPEN1, SEMOPEN2)
freqboth<-c(FREQ1, FREQ2)
allfreqboth<-c(ALLFREQ1, ALLFREQ2)
cor(semopenboth, freqboth, method="pearson") # -0.3462541  : medium strength negative correlation
plot(semopenboth, freqboth)
cor(semopenboth, allfreqboth, method="pearson") #-0.4026855
 #Version of model with non-relativized SEMOPEN
##C value  = 0.7187445 for main effects only model

 #Model with non-relativized SemOpen has HUGE colinearity with GENFREQDIFF 
 #Version of model with non-relativized logged SEMOPEN AND no GFdiff was worse than model with relativized SEMOPEN and GFdiff  (model2.glm6)
#C = 0.7880309 : lower than one w GenFreq and separate relativized SemOpen
#R2 = 0.351 not as good as other model, (model with GFdiff and rel SemOpen is 0.367)
#It did though have slightly higher percentage of correct predicted values  71% instead of 70%

 #instead, try Stefan's residuals trick to get the information from SemOpen that isn't due to GFdiff
SEMOPEN_nr<-c(mod1.ntype, mod2.ntype)
GENFREQ_nr<-c(mod1.freq, mod2.freq)
SeOpGF<-lm(SEMOPEN_nr ~GENFREQ_nr)
summary(SeOpGF)
hist(SeOpGF$residuals)
SemOpwoGF<-(SeOpGF$residuals)
SemOpwoGF<-unname(SemOpwoGF)
SemOpwoGF1<-c(SemOpwoGF[1:1810])
SemOpwoGF2<-c(SemOpwoGF[1811:3620])
SemOpwoGFdiff<-c(SemOpwoGF1-SemOpwoGF2, SemOpwoGF2-SemOpwoGF1)
SemOpwoGFdiff_orig<-c(SemOpwoGF1-SemOpwoGF2)
SemOpwoGFdiff_orig[which(SemOpdiff_orig>0.5)] #shows that some make more sense , like "black and red uniform" [1792] is not so huge of a difference now, "black and big bird" [1676] is still on the high end with 0.4183 (where hist shows most are btw -.25 and .25), but still it's below .5 now
#model with residuals has C = 0.7924258 :  better than other non-relativized SemOp and no GF, but lower than one w GenFreq and separate relativized SemOpen (0.7949879)
library(rms)
lrm(formula(model.glm7)) #R2 = 0.362 not as good as other model, (model with GFdiff and rel SemOpen is 0.367)
actual<-data.frame(RLdiff, SemOpwoGFdiff, ICdiff, SubObjdiff, ALdiff, NSpFdiff, GFdiff, CON)
classifications<-predict(model.glm7, type="response")
classif2<-ifelse(classifications>=0.5, "predict", "nonpredict")
evaluation<-table(classif2, ORDER)
addmargins(evaluation);sum(diag(evaluation))/sum(evaluation) #70.82873% correct
library(effects)
par(cex.main=0.1, cex.axis=0.5)
plot(allEffects(model.glm7), ask=FALSE, grid=TRUE, ylim=c(0,1), rescale.axis=FALSE, cex=1.0) 
#plotting effects separately
plot(effect("NSpFdiff:CON", model2.glm6))
plot(effect("SubObjdiff:CON", model2.glm6))


#GEN FREQ with ALL mod positions
#Get MOD Freq for MOD ROOTs in ANY lexeme/construction
# kaist.files <- dir(".", full.names=TRUE) #
# kaist.mod.matches <-c()
# for (i in 1:length(kaist.files)){
  # cat(i, "\n")
  # curr.kaist.file<-scan(kaist.files[i], what=character(0), sep="\n", quiet=TRUE)
  # # remove metadata lines
  # curr.kaist.file<-gsub("<.*", "", curr.kaist.file, perl=TRUE)
   # curr.mod.matches<-exact.matches("[^\\s]+\\s+[^\\s]+paa[^\\s]+", curr.kaist.file, pcre=TRUE, case.sens=FALSE)[[1]] #
  # # add to previous matches of adjective doublets 
  # kaist.mod.matches <- c(kaist.mod.matches, curr.mod.matches)
# } 
# length(kaist.mod.matches);head(kaist.mod.matches) #838445
# write.table(kaist.mod.matches, file="../KAIST_modall_orig.csv", quote=F, sep="\t", row.names=F, col.names=F)
kaist.mod.matches <-read.table(file=("../KAIST_modall_orig.csv"), quote="", sep="\t", header=F, comment.char="")
#map to kaist.mod (adj lexeme), and modanalysis (adj analysis to get root)
kaist.mod<-kaist.mod.matches$V1
kaist.modanalysis<-kaist.mod.matches$V2
#clean up weird punctuation and slang tags
kaist.modanalysis <-gsub(".\\/sl\\+", "", kaist.modanalysis, perl=TRUE) #clean up /sl (slang) tags
kaist.modanalysis <-gsub("\\+.\\/sr", "", kaist.modanalysis, perl=TRUE) #clean up /sl (slang) tags
kaist.modanalysis <-gsub("\\(.*\\)", "", kaist.modanalysis, perl=TRUE);length(kaist.modanalysis) #get rid of anything within parentheses in analysis columns
kaist.modanalysis <-gsub("[\\\\\\,&}\\.\\(\\)\\\"\\[\\]]+", "", kaist.modanalysis, perl=TRUE) #Get rid of rid of weird punctuation ["()[]\,&.]
#same thing for kaist.mod
kaist.mod <-gsub(".\\/sl\\+", "", kaist.mod, perl=TRUE) #clean up /sl (slang) tags
kaist.mod <-gsub("\\+.\\/sr", "", kaist.mod, perl=TRUE) #clean up /sl (slang) tags
kaist.mod <-gsub("\\(.*\\)", "", kaist.mod, perl=TRUE);length(kaist.mod) #get rid of anything within parentheses in analysis columns
kaist.mod <-gsub("[\\\\\\,&}\\.\\(\\)\\\"\\[\\]]+", "", kaist.mod, perl=TRUE) #Get rid of rid of weird punctuation ["()[]\,&.]
#get the adjective roots: 
kaist.modroot<-gsub("\\/[^\\s]+", "", kaist.modanalysis, perl=TRUE);head(kaist.modroot)
kaist.modroot<-gsub("[\\d\\w\\s]", "", kaist.modroot, perl=TRUE)

#Get GENFREQ for ALL positions for MOD1 and MOD2
	# # MOD1 Frequency
# mod1.allfreq<-c()
# for (j in 1:length(MOD1.ROOT)){ 
		# mod1.allfreq <-c(mod1.allfreq, length(which(kaist.modroot==MOD1.ROOT[j])))	
	# }
# #plus 1st position mod in mod mod n
# head(mod1.allfreq);length(mod1.allfreq)#1810
# unique(MOD1.ROOT[which(mod1.allfreq==(max(mod1.allfreq)))])
# hist(mod1.allfreq)
# plot(sort(table(mod1.allfreq)))
#write.table(mod1.allfreq, file="../KAIST_mod1.allfreq.csv", row.names=F, col.names=F,quote=F, sep="\t")
mod1.allfreq<-read.table(file="../KAIST_mod1.allfreq.csv", header=FALSE, sep="\t", quote="")
mod1.allfreq<-mod1.allfreq$V1
	
	# # MOD2 Frequency
mod2.allfreq<-c()
for (j in 1:length(MOD2.ROOT)){ 
		mod2.allfreq <-c(mod2.allfreq, length(which(kaist.modroot ==MOD2.ROOT[j])))	
	}
head(mod2.allfreq);length(mod2.allfreq)#1810
unique(MOD2.ROOT[which(mod2.allfreq==(max(mod2.allfreq)))])
hist(mod2.allfreq)
plot(sort(table(mod2.allfreq)))
mod2.allfreq[949]<-4 #fix this 차이나식 one
write.table(mod2.allfreq, file="../KAIST_mod2.allfreq.csv", quote=F, row.names=F, col.names=F, sep="\t")
#mod2.freq<-read.table(file="../KAIST_mod2.freq.csv", header=FALSE, sep="\t", quote="")
#mod2.freq<-mod2.freq$V1

#Calculate Semantic openness
SEMOPEN_both<-c(mod1.ntype/mod1.allfreq, mod2.ntype/mod2.allfreq)
# # #Prediction SEMOPEN1 > SEMOPEN2 NOUN
SemOpdiff <-c((mod1.ntype/mod1.allfreq)-(mod2.ntype/mod2.allfreq), (mod2.ntype/mod2.allfreq)-(mod1.ntype/mod1.allfreq))
##SEM OPEN (END)

#GENFREQ (START)
GFboth<-c(mod1.allfreq, mod2.allfreq)
GFdiff<-c((log(mod1.allfreq)-log(mod2.allfreq)),(log(mod2.allfreq)-log(mod1.allfreq)))
#GENFREQ (END)
#reconstructing mod1.comp and mod2.comp
	#INDCOMP1<-mod1.comp/mod1.freq
	#INDCOMP2<-mod2.comp/mod2.freq
	#so multiplying INDCOMP1 by mod1.freq will result in mod1.comp
mod1.comp<-c(INDCOMP1*mod1.freq)
mod2.comp<-c(INDCOMP2*mod2.freq)
ICdiff<-c((mod1.comp/mod1.allfreq)-(mod2.comp/mod2.allfreq), (mod2.comp/mod2.allfreq)-(mod1.comp/mod1.allfreq))

#Updating kaist.everything table with new values for GF, SEMOPEN and INDCOMP
# SEMOPEN1<-c(mod1.ntype/mod1.allfreq)
# SEMOPEN2<-c(mod2.ntype/mod2.allfreq)
# ALLFREQ1<-c(mod1.allfreq)
# ALLFREQ2<-c(mod2.allfreq)
# INDCOMP1<-c(mod1.comp/mod1.allfreq)
# INDCOMP2<-c(mod2.comp/mod2.allfreq)
#kaist.everything2<-data.frame(FILE, PHRASE, PHRASE_DEF, MOD1, MOD1.ROOT, MOD1_ANALYSIS, MOD1_DEF, CONNECT, CONNECT2, RT.LENGTH1, SEMOPEN1, INDCOMP1, SUBJOBJ1, MOD1_AFFLOAD, NSPEC1, NTYPE1, COMP1, ALLFREQ1, FREQ1, MOD2, MOD2.ROOT, MOD2_ANALYSIS, MOD2_DEF, RT.LENGTH2, SEMOPEN2, INDCOMP2, SUBJOBJ2, MOD2_AFFLOAD, NSPEC2, NTYPE2, COMP2, ALLFREQ2, FREQ2, NOUN, NOUN.ROOT, NOUN_ANALYSIS, NOUNFREQ)
#write.table(kaist.everything2, file="../KAIST_everything2.csv", quote=F, sep="\t", row.names=F, col.names=T)
kaisteverything2<-read.table(file="/Users/heathersimpson/Documents/Corpora/KAIST/KAIST_everything2.csv", quote="", sep="\t", header=T)
attach(kaisteverything2)


##REGRESSION ANALYSIS
RTLENGTHDIFF<-c(RT.LENGTH2-RT.LENGTH1, RT.LENGTH1-RT.LENGTH2);length(RTLENGTHDIFF)
SEMOPENDIFF<-c(SEMOPEN1-SEMOPEN2, SEMOPEN2-SEMOPEN1)
INDCOMPDIFF<-c(INDCOMP1-INDCOMP2, INDCOMP2-INDCOMP1)
SUBJOBJDIFF<-c(SUBJOBJ1-SUBJOBJ2, SUBJOBJ2-SUBJOBJ1)
AFFLOADDIFF<-c(MOD1_AFFLOAD-MOD2_AFFLOAD, MOD2_AFFLOAD-MOD1_AFFLOAD)
NSPECDIFF<-c(NSPEC2-NSPEC1,NSPEC1-NSPEC2)
GENFREQDIFF<-c(log(ALLFREQ1)-log(ALLFREQ2),log(ALLFREQ2)-log(ALLFREQ1))
dep<-c(rep("predict", 1810), rep("nonpredict", 1810))
dep<-as.factor(dep)
#Variable for connecting grammar
CONNECTDIFF<-append(CONNECT2, CONNECT2)
CONNECTDIFF <-factor(CONNECTDIFF, levels=1:nlevels(CONNECT2), labels=levels(CONNECT2))
table(CONNECTDIFF)


#Make version with shorter names for plot
RLdiff<-RTLENGTHDIFF
RLdiff<-as.vector(RLdiff)
SemOpdiff<-SEMOPENDIFF
SemOpdiff<-as.vector(SemOpdiff)
ICdiff<-INDCOMPDIFF
ICdiff<-as.vector(ICdiff)
SubObjdiff<-SUBJOBJDIFF
SubObjdiff<-as.vector(SubObjdiff)
ALdiff<-AFFLOADDIFF
ALdiff<-as.vector(ALdiff)
NSpFdiff<-NSPECDIFF
NSpFdiff<-as.vector(NSpFdiff)
GFdiff<-GENFREQDIFF
GFdiff<-as.vector(GFdiff)
CON<-CONNECTDIFF
ORDER<-dep

#Change CON levels: 
CON2<- factor(CON, levels = c(levels(CON), "noncoord", "coord"))
CON2[CON2=="etm"]<-"noncoord"
CON2[CON2=="ecc"]<-"coord"

#Generalized Linear Regression Model Analysis with ALLFREQ measures
options(contrasts=c("contr.treatment", "contr.poly"))
model.glm<-glm(ORDER ~ (RLdiff+ SemOpdiff + ICdiff + SubObjdiff + ALdiff + NSpFdiff +  GFdiff + CON), family=binomial)
library(Hmisc)
somers2(binomial()$linkinv(fitted(model.glm)), as.numeric(dep)-1) # C: 0.7201438
library(car)
vif(model.glm) # higest is GFdiff with 2.129766
add1(model.glm, scope= ~((RLdiff+ SemOpdiff + ICdiff + SubObjdiff + ALdiff + NSpFdiff + GFdiff +CON)^2), test=c("Chisq"))
model.glm2<-glm(ORDER ~ (RLdiff+ SemOpdiff + ICdiff + SubObjdiff + ALdiff + NSpFdiff +  GFdiff + CON+RLdiff:CON), family=binomial) #add RLdiff:CON
add1(model.glm2, scope= ~((RLdiff+ SemOpdiff + ICdiff + SubObjdiff + ALdiff + NSpFdiff + GFdiff +CON)^2), test=c("Chisq")) 
model.glm3<-glm(ORDER ~ (RLdiff+ SemOpdiff + ICdiff + SubObjdiff + ALdiff + NSpFdiff +  GFdiff + CON+RLdiff:CON+ SubObjdiff:CON), family=binomial) #SubObjdiff:CON
add1(model.glm3, scope= ~((RLdiff+ SemOpdiff + ICdiff + SubObjdiff + ALdiff + NSpFdiff + GFdiff +CON)^2), test=c("Chisq")) 
model.glm4<-glm(ORDER ~ (RLdiff+ SemOpdiff + ICdiff + SubObjdiff + ALdiff + NSpFdiff +  GFdiff + CON+RLdiff:CON+ SubObjdiff:CON + NSpFdiff:CON), family=binomial) #NSpFdiff:CON
add1(model.glm4, scope= ~((RLdiff+ SemOpdiff + ICdiff + SubObjdiff + ALdiff + NSpFdiff + GFdiff +CON)^2), test=c("Chisq")) 
model.glm5<-glm(ORDER ~ (RLdiff+ SemOpdiff + ICdiff + SubObjdiff + ALdiff + NSpFdiff +  GFdiff + CON+RLdiff:CON+ SubObjdiff:CON + NSpFdiff:CON+ ICdiff:CON), family=binomial) #ICdiff:CON
add1(model.glm5, scope= ~((RLdiff+ SemOpdiff + ICdiff + SubObjdiff + ALdiff + NSpFdiff + GFdiff +CON)^2), test=c("Chisq")) 
model.glm6<-glm(ORDER ~ (RLdiff+ SemOpdiff + ICdiff + SubObjdiff + ALdiff + NSpFdiff +  GFdiff + CON+RLdiff:CON+ SubObjdiff:CON + NSpFdiff:CON+ ICdiff:CON+ SemOpdiff:CON), family=binomial) #SemOpdiff:CON
add1(model.glm6, scope= ~((RLdiff+ SemOpdiff + ICdiff + SubObjdiff + ALdiff + NSpFdiff + GFdiff +CON)^2), test=c("Chisq")) 
model.glm7<-glm(ORDER ~ (RLdiff+ SemOpdiff + ICdiff + SubObjdiff + ALdiff + NSpFdiff +  GFdiff + CON+RLdiff:CON+ SubObjdiff:CON + NSpFdiff:CON+ ICdiff:CON+ SemOpdiff:CON + GFdiff:CON), family=binomial) #GFdiff:CON -significant p-value and lowers AIC
add1(model.glm7, scope= ~((RLdiff+ SemOpdiff + ICdiff + SubObjdiff + ALdiff + NSpFdiff + GFdiff +CON)^2), test=c("Chisq")) 
summary(model.glm7)
# Deviance Residuals: 
    # Min       1Q   Median       3Q      Max  
# -2.7042  -0.9554   0.0000   0.9554   2.7042  

# Coefficients:
                    # Estimate Std. Error z value Pr(>|z|)    
# (Intercept)       -2.513e-15  4.104e-02   0.000   1.0000    
# RLdiff             5.211e-01  5.251e-02   9.923  < 2e-16 ***
# SemOpdiff         -1.821e+00  3.620e-01  -5.030 4.91e-07 ***
# ICdiff            -1.998e+00  2.400e+00  -0.833   0.4050    
# SubObjdiff        -1.011e-01  1.664e-02  -6.079 1.21e-09 ***
# ALdiff             2.708e-01  6.797e-02   3.984 6.78e-05 ***
# NSpFdiff           1.993e+01  1.896e+00  10.510  < 2e-16 ***
# GFdiff             1.102e-01  2.709e-02   4.068 4.75e-05 ***
# CONetm             1.467e-16  1.278e-01   0.000   1.0000    
# RLdiff:CONetm     -1.185e+00  1.535e-01  -7.717 1.19e-14 ***
# SubObjdiff:CONetm  1.803e-01  4.181e-02   4.313 1.61e-05 ***
# NSpFdiff:CONetm    7.134e+01  1.257e+01   5.674 1.40e-08 ***
# ICdiff:CONetm      6.229e+01  8.042e+00   7.746 9.49e-15 ***
# SemOpdiff:CONetm   4.371e+00  9.672e-01   4.519 6.21e-06 ***
# GFdiff:CONetm     -1.766e-01  8.452e-02  -2.089   0.0367 *  
library(car)
options(contrasts=c("contr.sum", "contr.poly"))
Anova(model.glm7, type="III", test.statistic="Wald")# results in just about exactly the same as above
library(rms)
lrm(formula(model.glm7)) #R2 = 0.364, C= 0.797

#final model without SEMOPEN
#model.glm6<-glm(ORDER ~ (RLdiff + ICdiff + SubObjdiff + ALdiff + NSpFdiff +  GFdiff + CON+RLdiff:CON+ ICdiff:CON + NSpFdiff:CON+ SubObjdiff:CON+ GFdiff:CON), family=binomial) #R2 = 0.355, C = .791

actual<-data.frame(RLdiff, SemOpdiff, ICdiff, SubObjdiff, ALdiff, NSpFdiff, GFdiff, CON)
classifications<-predict(model.glm7, type="response")
classif2<-ifelse(classifications>=0.5, "predict", "nonpredict")
evaluation<-table(classif2, ORDER)
addmargins(evaluation);sum(diag(evaluation))/sum(evaluation) #71.5% correct
library(effects)
par(cex.main=0.1, cex.axis=0.5)
plot(allEffects(model.glm7), ask=FALSE, grid=TRUE, ylim=c(0,1), rescale.axis=FALSE, cex=1.0) 
#plotting effects separately
plot(effect("RLdiff:CON", model.glm7))
plot(effect("ICdiff:CON", model.glm7))
plot(effect("SemOpdiff:CON", model.glm7))
plot(effect("GFdiff:CON", model.glm7))
plot(effect("NSpFdiff:CON", model.glm7))
plot(effect("SubObjdiff:CON", model.glm7))
plot(effect("RLdiff", model.glm7))
plot(effect("SubObjdiff", model.glm7))

#ncps vs. paa
kaistncps2<-read.table(file="../KAIST_ncps2.csv", quote="", sep="\t", header=F)
allmod<-c(as.character(MOD1), as.character(MOD2))
allmodanalysis<-c(as.character(MOD1_ANALYSIS), as.character(MOD2_ANALYSIS))
allpaahada<-unique(allmod[grep("하\\/", allmodanalysis)])
allncps <-unique(kaistncps2$V1[grep("ncps", kaistncps2$V2)])
length(allpaahada) #113
length(allncps) #386
write.table(allpaahada, file="../KAIST_allpaahada.txt", quote=F, sep="\t", row.names=F, col.names=F)
write.table(allncps, file="../KAIST_allncps.txt", quote=F, sep="\t", row.names=F, col.names=F)

#TRYING TO FIGURE OUT WHAT PERCENTAGE OF EACH ADJECTIVE'S NOUN COMBOS ARE HAPAXES
#kaist.modn.modroot
#kaist.modn.roots
# length(unique(kaist.modn.modroot)) #1830
# length(unique(kaist.modn.roots)) #117202
# x<-hist(table(kaist.modn.roots)) #719 files appeared between 0 and 5 times
# 719/780 #92.2%
##calculating modn.freq = for each modn root combo (non unique), how many items in the vector of modn root combos are equal to it, i.e. its corpus frequency, = length of kaist.modn.modroots and kaist.modn.roots
##if I was going to be totally accurate I need to add mod1 combos to kaist.modn.roots also
# modn.freq <-c()
# for (j in 1:length(kaist.modn.roots)){ 
		# modn.freq <-c(modn.freq, length(which(kaist.modn.roots==kaist.modn.roots[j])))	
	# }
 head(modn.freq);length(modn.freq)#
 #write.table(modn.freq, file="../KAIST_modnfreq.csv", quote=F, sep="\t", row.names=F, col.names=F)
modn.freq<-read.table(file="../KAIST_modnfreq.csv", quote="", sep="\t", header=T)
modn.freq<-c(modn.freq$x)
 unique.modroots<-as.character(unique(as.character(MOD1.ROOT, MOD2.ROOT))) #length(unique.modroots) = 234
 mod.roots.both<-c(as.character(MOD1.ROOT), as.character(MOD2.ROOT))

 unique.modn.counts<-unique(data.frame(kaist.modn.modroot, kaist.modn.roots, modn.freq))
 unique.kaist.mod.roots<-as.character(unique.modn.counts$kaist.modn.modroot)
 unique.kaist.modn.roots<-as.character(unique.modn.counts$kaist.modn.roots)
 unique.modn.freq<-c(unique.modn.counts$modn.freq)
 rm(unique.modn.counts)
 modncombo.freq<-list()
for (j in 1:length(mod.roots.both)){ 
		modncombo.freq[j]<-list(c(unique.modn.freq[which(unique.kaist.mod.roots == mod.roots.both[j])]))	
	}
length(modncombo.freq);head(modncombo.freq)
length(modncombo.freq[[6]]==1)
#calculate for each list position, the percentage of adj+noun collocation types that only occured once
hapaxperc<-c()
for(i in 1:length(modncombo.freq)){
	hapaxperc<-c(hapaxperc, sum(modncombo.freq[[i]]==1)/length(modncombo.freq[[i]]))
} #at first I did sum(modncombo.freq[[i]]), but that gets the percentage of adj+noun collocation tokens, I want percentage of types
 allfreq.both<-c(ALLFREQ1, ALLFREQ2)
cor(hapaxperc, log(allfreq.both), use="complete.obs", method=c("pearson")) #r = -0.576 (unlogged allfreq.both is -0.364)
#use "complete.obs" to remove 3 NaN values (created when the # of matches in unique.kaist.mod.roots == 0 , this happened for mod.roots.both[c(329, 918, 1053)], "어줍", "허옇", "거쿨지", these must have only occured in mod1 attr. position? 
plot(hapaxperc, log(allfreq.both))
#Now try with attributive only frequency
 freq.both<-c(FREQ1, FREQ2)
cor(hapaxperc, log(freq.both), use="complete.obs", method=c("pearson"))  #r = -0.6552342
#Now try with unique mods : 
hapallf<-unique(data.frame(mod.roots.both, hapaxperc, allfreq.both))
str(hapallf) #352 obs for 352 unique mods
cor(hapallf$hapaxperc, log(hapallf$allfreq.both), use="complete.obs", method=c("pearson")) #r = -0.3956076

#Investigating SUBJOBJ, getting Examples for paper
#checking on SUBJOBJDIFF
PHRASE_DEF[(which(SUBJOBJ1-SUBJOBJ2<0))] 
SUBJOBJcombos<-paste(SUBJOBJ1, SUBJOBJ2, sep=" ")
subjobj_diffonly<-
sort(table(SUBJOBJcombos))
#most common of ones that have a difference: 
#11 12 with 115
#11 5 with 94
#5 11 with 61  
#note 11 11, 12 12, 8 8, 5 5  are most common of equal ones
#12 11 with 48 
# 11 10 with 41  
# 7 12  with 40, 						12 7 has 32 
#young beautiful woman 
#picking out some ones that sound terrible to me as English speaker
#black and soft moss
#straight and long fur
#red and long mustache
#light and cool-looking bird 
#small and cute pet
#green and big plant
#cold and big room
#cheap and good stuff
#hot and long summer
#Could "and" be blocking the SUBJOBJ? 
#Proportion of CONNECT with SUBJOBJ<0
#investgating patterns in SUBJOBJ diff and CONNECT
#are there any combos that are in "etm" set that are not in "ecc" set? 
setdiff(SUBJOBJcombos[CONNECT2=="etm"],SUBJOBJcombos[CONNECT2=="ecc"])
#answer: YES, which could be strong evidence for lexical preference effect of some type since etm is much less frequent than ecc
	#"12 9"  "8 7"   "9 7"   "7 5"   "10 13" "8 5"
#"12 9" ,  only one example: "strong and fast songs"
#  "8 7"  
	# [1] strange and young knight      sad and young couple          Non-sensical and young person
	# [4] mischievous and young person  sly and young oriental-doctor sly and young oriental-doctor
	# [7] unique and new fact           happy and young child         energetic and young man-power
	#[10] energetic and young man-power mean and young guy            happy and young moment       
	#[13] sad and young moment 
#"9 7"   
	#[1] fast and new boat-road
#"7 5"
	#[1] old and black-ish tree
#"10 13"
	#[1] damp and unfamiliar place
#"8 5"
	#[1] expensive and blue letter-sheet
par(mfrow=c(2, 1))
plot(prop.table(table(SUBJOBJcombos[CONNECT2=="ecc"])), ylab=c("ecc"))
plot(prop.table(table(SUBJOBJcombos[CONNECT2=="etm"])), ylab=c("etm"))
par(mfrow=c(1, 1))
#11 5 , 6 11, and 8 11 seems much higher percentage for etm than for ecc 
#How about just do a Chi square test for whether these two affect each other 
SUBJOBJcombos2<-as.factor(SUBJOBJcombos)
SubObjCon<-table(SUBJOBJcombos2, CONNECT2)
eval.SubObjCon<-chisq.test(SubObjCon, correct=F); eval.SubObjCon
assocplot(SubObjCon)
x<-assocplot(SubObjCon)
#expected frequencies
eval.SubObjCon$expected
#computing Cramer's V effect size
sqrt((eval.SubObjCon$statistic)/sum(eval.SubObjCon)*(min(dim(eval.SubObjCon))-1)) # doesn't work, gives error message about list type
#compute contribution of variable level combinations
eval.SubObjCon$residuals^2

#what if I just do chisq for the levels of SUBJOBJ for each connect type separately? 
soecc1<-as.factor(SUBJOBJ1[CONNECT2=="ecc"])
soecc2<-as.factor(SUBJOBJ2[CONNECT2=="ecc"])
SubObjecc<-table(soecc1, soecc2)
eval.SubObjecc<-chisq.test(SubObjecc, correct=F); eval.SubObjecc
eval.SubObjecc$expected
eval.SubObjecc$residuals^2

soetm1<-as.factor(SUBJOBJ1[CONNECT2=="etm"])
soetm2<-as.factor(SUBJOBJ2[CONNECT2=="etm"])
SubObjetm<-table(soetm1, soetm2)
eval.SubObjetm<-chisq.test(SubObjetm, correct=F); eval.SubObjetm
eval.SubObjetm$expected
eval.SubObjetm$residuals^2

#There are so many levels that it might be problematic with the expected frequencies, maybe try just for the modifier the connector is on
SUBJOBJ1_2<-as.factor(SUBJOBJ1)
SubObjCon<-table(SUBJOBJ1_2, CONNECT2)
eval.SubObjCon<-chisq.test(SubObjCon, correct=F); eval.SubObjCon
sqrt((eval.SubObjCon$statistic)/sum(eval.SubObjCon)*(min(dim(eval.SubObjCon))-1)) # doesn't work, gives error message about list type
#compute contribution of variable level combinations
eval.SubObjCon$residuals^2
          # CONNECT2
# SUBJOBJ1_2          etm          ecc
        # 5   4.352192608  0.868131196
        # 6  10.725758419  2.139465397
        # 7   1.671152936  0.333344622
        # 8   3.734034083  0.744827209
        # 9   0.895361074  0.178597537
        # 10  0.091367966  0.018225154
        # 11  2.060397739  0.410987223
        # 12  3.058392342  0.610057054
        # 13  0.023129689  0.004613675
assocplot(SubObjCon)
#more (6, 8, and 11) in etm, less of (5, 7, 9, 10, 12, 13)

adjpos<-c(rep("1", 1810), rep("2", 1810))
subobjpos<-table(subobjboth, adjpos)
eval.subobjpos<-chisq.test(subobjpos, correct=F); eval.subobjpos #p-value = 2.337e-07
eval.subobjpos$residuals^2

#GET EXAMPLES OF EACH VARIABLE FOR PAPER
#get versions of variable with both mod1 and mod2 values included (to get more accurate distribution)
rlboth<-c(RT.LENGTH1, RT.LENGTH2)
semopboth<-c(SEMOPEN1, SEMOPEN2)
icboth<-c(INDCOMP1, INDCOMP2)
subobjboth<-c(SUBJOBJ1, SUBJOBJ2)
alboth<-c(MOD1_AFFLOAD, MOD2_AFFLOAD)
AFFLOAD1<-MOD1_AFFLOAD
AFFLOAD2<-MOD2_AFFLOAD
nspfboth<-c(NSPEC1, NSPEC2)
gfboth<-c(ALLFREQ1, ALLFREQ2)
mod.roots.both<-c(as.character(MOD1.ROOT), as.character(MOD2.ROOT))
mod.def.both<-c(as.character(MOD1_DEF), as.character(MOD2_DEF))
modanalysis.both<-c(as.character(MOD1_ANALYSIS), as.character(MOD2_ANALYSIS))
nounanalysis.both<-c(as.character(NOUN_ANALYSIS), as.character(NOUN_ANALYSIS))
phrase.both<-c(as.character(PHRASE), as.character(PHRASE))
phrase.def.both<-c(as.character(PHRASE_DEF), as.character(PHRASE_DEF))
compboth<-c(COMP1, COMP2)
nfreqboth<-c(   cv, NOUNFREQ)

#INDCOMP examples
hist(icboth)
boxplot(icboth)
plot(sort(icboth))
plot(sort(log(icboth)), log(1:length(icboth)))

sum(icboth<=0.01)
(mod.roots.both[icboth<=0.01])[222]
(mod.def.both[icboth<=0.01])[222]
(gfboth[icboth<=0.01])[222]
(compboth[icboth<=0.01])[222]
(icboth[icboth<=0.01])[222]

sum(icboth<0.01 & icboth>0.005)
(mod.roots.both[icboth<0.01 & icboth>0.005])[222]
(mod.def.both[icboth<0.01 & icboth>0.005])[222]
(gfboth[icboth<0.01 & icboth>0.005])[222]
(compboth[icboth<0.01 & icboth>0.005])[222]
(icboth[icboth<0.01 & icboth>0.005])[222]

sum(icboth<0.02 & icboth>0.01)
(mod.roots.both[icboth<0.02 & icboth>0.01])[222]
(mod.def.both[icboth<0.02 & icboth>0.01])[222]
(gfboth[icboth<0.02 & icboth>0.01])[222]
(compboth[icboth<0.02 & icboth>0.01])[222]
(icboth[icboth<0.02 & icboth>0.01])[222]

sum(icboth<0.03 & icboth>0.02)
(mod.roots.both[icboth<0.03 & icboth>0.02])[222]
(mod.def.both[icboth<0.03 & icboth>0.02])[222]
(gfboth[icboth<0.03 & icboth>0.02])[222]
(compboth[icboth<0.03 & icboth>0.02])[222]
(icboth[icboth<0.03 & icboth>0.02])[222]

sum(icboth<0.05 & icboth>0.03)
(mod.roots.both[icboth<0.05 & icboth>0.03])[222]
(mod.def.both[icboth<0.05 & icboth>0.03])[222]
(gfboth[icboth<0.05 & icboth>0.03])[222]
(compboth[icboth<0.05 & icboth>0.03])[222]
(icboth[icboth<0.05 & icboth>0.03])[222]

sum(icboth>=0.06)
(mod.roots.both[icboth>=0.06])[28]
(mod.def.both[icboth>=0.06])[28]
(gfboth[icboth>=0.06])[28]
(compboth[icboth>=0.06])[28]
(icboth[icboth>=0.06])[28]

#NSPECFREQ examples
modncolfreqboth<-c(mod1nmatches2, mod2nmatches2) 
modnboth<-c(paste(MOD1, NOUN, sep=" "), paste(MOD2, NOUN, sep=" "))
# def.split<-strsplit(as.vector(PHRASE_DEF), "\\s", perl=TRUE)
# NOUN_DEF<-sapply(def.split, "[", 3);head(NOUN_DEF)
# modndefboth<-c(paste(MOD1_DEF, NOUN, sep=" "), paste(MOD1_DEF, NOUN, sep=" "))

hist(nspfboth)
boxplot(nspfboth)
plot(sort(nspfboth))
plot(sort(log(nspfboth)), log(1:length(nspfboth)))
x<-hist(nspfboth);x

sum(nspfboth<=0.05) #3369
(modnboth[nspfboth<= 0.05])[222]
(phrase.def.both[nspfboth<= 0.05])[222]
(modncolfreqboth[nspfboth<= 0.05])[222]
(nfreqboth[nspfboth<= 0.05])[222]
(nspfboth[nspfboth<= 0.05])[222]

sum(nspfboth>0.05 & nspfboth<0.10) #100
(modnboth[nspfboth>0.05 & nspfboth<0.10])[22]
(phrase.def.both[nspfboth>0.05 & nspfboth<0.10])[22]
(modncolfreqboth[nspfboth>0.05 & nspfboth<0.10])[22]
(nfreqboth[nspfboth>0.05 & nspfboth<0.10])[22]
(nspfboth[nspfboth>0.05 & nspfboth<0.10])[22]

sum(nspfboth>0.10 & nspfboth<0.2)#49
(modnboth[nspfboth>0.10 & nspfboth<0.2])[22]
(phrase.def.both[nspfboth>0.10 & nspfboth<0.2])[22]
(modncolfreqboth[nspfboth>0.10 & nspfboth<0.2])[22]
(nfreqboth[nspfboth>0.10 & nspfboth<0.2])[22]
(nspfboth[nspfboth>0.10 & nspfboth<0.2])[22]

sum(nspfboth>0.2 & nspfboth<0.4)#34
(modnboth[nspfboth>0.2 & nspfboth<0.4])[18]
(phrase.def.both[nspfboth>0.2 & nspfboth<0.4])[18]
(modncolfreqboth[nspfboth>0.2 & nspfboth<0.4])[18]
(nfreqboth[nspfboth>0.2 & nspfboth<0.4])[18]
(nspfboth[nspfboth>0.2 & nspfboth<0.4])[18]

sum(nspfboth>0.4 & nspfboth<0.5) #2
(modnboth[nspfboth>0.4 & nspfboth<0.5])
(phrase.def.both[nspfboth>0.4 & nspfboth<0.5])
(modncolfreqboth[nspfboth>0.4 & nspfboth<0.5])
(nfreqboth[nspfboth>0.4 & nspfboth<0.5])
(nspfboth[nspfboth>0.4 & nspfboth<0.5])

sum(nspfboth>=0.6) #36
(modnboth[nspfboth>=0.6])[28]
(phrase.def.both[nspfboth>=0.6])[28]
(modncolfreqboth[nspfboth>=0.6])[28]
(nfreqboth[nspfboth>=0.6])[28]
(nspfboth[nspfboth>=0.6])[28]


#GENFREQ examples
hist(gfboth)
boxplot(gfboth)

#GET EXAMPLES OF ORDERINGS FOR PAPER
#Get versions of each variable with only PREDICT values
RLdiff_orig<-c(RT.LENGTH2-RT.LENGTH1)
SemOpdiff_orig<-c(SEMOPEN1-SEMOPEN2)
ICdiff_orig<-c(INDCOMP1-INDCOMP2)
SubObjdiff_orig<-c(SUBJOBJ1-SUBJOBJ2)
ALdiff_orig<-c(MOD1_AFFLOAD-MOD2_AFFLOAD)
NSpFdiff_orig<-c(NSPEC2-NSPEC1)
GFdiff_orig<-c(log(ALLFREQ1)-log(ALLFREQ2))
CON_orig<-CONNECT2

 #INDCOMP diff ordering
hist(ICdiff_orig)
x<-hist(ICdiff_orig);x
 #expected

#unexpected
sum(ICdiff_orig>=-0.02 & ICdiff_orig<=0.00) #737
(PHRASE[ICdiff_orig>=-0.02 & ICdiff_orig<=0.00])[111]
(PHRASE_DEF[ICdiff_orig>=-0.02 & ICdiff_orig<=0.00])[111]
(INDCOMP1[ICdiff_orig>=-0.02 & ICdiff_orig<=0.00])[111]
(INDCOMP2[ICdiff_orig>=-0.02 & ICdiff_orig<=0.00])[111]

sum(ICdiff_orig>=-0.04 & ICdiff_orig<=-0.03) #35
(PHRASE[ICdiff_orig>=-0.04 & ICdiff_orig<=-0.03])[22]
(PHRASE_DEF[ICdiff_orig>=-0.04 & ICdiff_orig<=-0.03])[22]
(INDCOMP1[ICdiff_orig>=-0.04 & ICdiff_orig<=-0.03])[22]
(INDCOMP2[ICdiff_orig>=-0.04 & ICdiff_orig<=-0.03])[22]

sum(ICdiff_orig>=-0.08 & ICdiff_orig<=-0.06) #3
(PHRASE[ICdiff_orig>=-0.08 & ICdiff_orig<=-0.06])[2]
(PHRASE_DEF[ICdiff_orig>=-0.08 & ICdiff_orig<=-0.06])[2]
(INDCOMP1[ICdiff_orig>=-0.08 & ICdiff_orig<=-0.06])[2]
(INDCOMP2[ICdiff_orig>=-0.08 & ICdiff_orig<=-0.06])[2]

#AFFLOAD diff ordering
hist(ALdiff_orig)
x<-hist(ALdiff_orig);x
 #expected
sum(ALdiff_orig>0 & ALdiff_orig<=1) #162
(PHRASE[ALdiff_orig>0 & ALdiff_orig<=1])[111]
(PHRASE_DEF[ALdiff_orig>0 & ALdiff_orig<=1])[111]
(AFFLOAD1[ALdiff_orig>0 & ALdiff_orig<=1])[111]
(AFFLOAD2[ALdiff_orig>0 & ALdiff_orig<=1])[111]

sum(ALdiff_orig>0 & ALdiff_orig<=1) #162
(PHRASE[ALdiff_orig>0 & ALdiff_orig<=1])[22]
(PHRASE_DEF[ALdiff_orig>0 & ALdiff_orig<=1])[22]
(AFFLOAD1[ALdiff_orig>0 & ALdiff_orig<=1])[22]
(AFFLOAD2[ALdiff_orig>0 & ALdiff_orig<=1])[22]

sum(ALdiff_orig>1 & ALdiff_orig<=2) #67
(PHRASE[ALdiff_orig>1 & ALdiff_orig<=2])[24]
(PHRASE_DEF[ALdiff_orig>1 & ALdiff_orig<=2])[c(1:25)]
(AFFLOAD1[ALdiff_orig>1 & ALdiff_orig<=2])[c(1:25)]
(AFFLOAD2[ALdiff_orig>1 & ALdiff_orig<=2])[c(1:25)]
#unexpected
sum(ALdiff_orig>=-1 & ALdiff_orig<0) #196
(PHRASE[ALdiff_orig>=-1 & ALdiff_orig<0])[111]
(PHRASE_DEF[ALdiff_orig>=-1 & ALdiff_orig<0])[111]
(AFFLOAD1[ALdiff_orig>=-1 & ALdiff_orig<0])[111]
(AFFLOAD2[ALdiff_orig>=-1 & ALdiff_orig<0])[111]

sum(ALdiff_orig>=-1 & ALdiff_orig<0) #196
(PHRASE[ALdiff_orig>=-1 & ALdiff_orig<0])[22]
(PHRASE_DEF[ALdiff_orig>=-1 & ALdiff_orig<0])[22]
(AFFLOAD1[ALdiff_orig>=-1 & ALdiff_orig<0])[22]
(AFFLOAD2[ALdiff_orig>=-1 & ALdiff_orig<0])[22]

sum(ALdiff_orig>=-2 & ALdiff_orig<=-1.1) #9
(PHRASE[ALdiff_orig>=-2 & ALdiff_orig<=-1.1])[2]
(PHRASE_DEF[ALdiff_orig>=-2 & ALdiff_orig<=-1.1])[2]
(AFFLOAD1[ALdiff_orig>=-2 & ALdiff_orig<=-1.1])[2]
(AFFLOAD2[ALdiff_orig>=-2 & ALdiff_orig<=-1.1])[2]

 #NSPEC diff ordering
hist(NSpFdiff_orig)
x<-hist(NSpFdiff_orig);x
 #expected
sum(NSpFdiff_orig>=0 & NSpFdiff_orig<=0.1) #1107
(PHRASE[NSpFdiff_orig>=0 & NSpFdiff_orig<=0.1])[111]
(PHRASE_DEF[NSpFdiff_orig>=0 & NSpFdiff_orig<=0.1])[111]
(NSPEC1[NSpFdiff_orig>=0 & NSpFdiff_orig<=0.1])[111]
(NSPEC2[NSpFdiff_orig>=0 & NSpFdiff_orig<=0.1])[111]

sum(NSpFdiff_orig>=0.10 & NSpFdiff_orig<=0.3) #61
(PHRASE[NSpFdiff_orig>=0.10 & NSpFdiff_orig<=0.3])[58]
(PHRASE_DEF[NSpFdiff_orig>=0.10 & NSpFdiff_orig<=0.3])[58]
(NSPEC1[NSpFdiff_orig>=0.10 & NSpFdiff_orig<=0.3])[58]
(NSPEC2[NSpFdiff_orig>=0.10 & NSpFdiff_orig<=0.3])[58]

sum(NSpFdiff_orig>=0.90 & NSpFdiff_orig<=1.0) #33
(PHRASE[NSpFdiff_orig>=0.90 & NSpFdiff_orig<=1.0])[2]
(PHRASE_DEF[NSpFdiff_orig>=0.90 & NSpFdiff_orig<=1.0])[2]
(NSPEC1[NSpFdiff_orig>=0.90 & NSpFdiff_orig<=1.0])[2]
(NSPEC2[NSpFdiff_orig>=0.90 & NSpFdiff_orig<=1.0])[2]

#unexpected
sum(NSpFdiff_orig>=-0.05 & NSpFdiff_orig<=0.00) #629
(PHRASE[NSpFdiff_orig>=-0.05 & NSpFdiff_orig<=0.00])[500]
(PHRASE_DEF[NSpFdiff_orig>=-0.05 & NSpFdiff_orig<=0.00])[500]
(NSPEC1[NSpFdiff_orig>=-0.05 & NSpFdiff_orig<=0.00])[500]
(NSPEC2[NSpFdiff_orig>=-0.05 & NSpFdiff_orig<=0.00])[500]

sum(NSpFdiff_orig>=-0.1 & NSpFdiff_orig<=-0.05) #35
(PHRASE[NSpFdiff_orig>=-0.1 & NSpFdiff_orig<=-0.05])[22]
(PHRASE_DEF[NSpFdiff_orig>=-0.1 & NSpFdiff_orig<=-0.05])[22]
(NSPEC1[NSpFdiff_orig>=-0.1 & NSpFdiff_orig<=-0.05])[22]
(NSPEC2[NSpFdiff_orig>=-0.1 & NSpFdiff_orig<=-0.05])[22]

sum(NSpFdiff_orig>=-0.3 & NSpFdiff_orig<=-0.1) #9
(PHRASE[NSpFdiff_orig>=-0.3 & NSpFdiff_orig<=-0.1])[2]
(PHRASE_DEF[NSpFdiff_orig>=-0.3 & NSpFdiff_orig<=-0.1])[2]
(NSPEC1[NSpFdiff_orig>=-0.3 & NSpFdiff_orig<=-0.1])[2]
(NSPEC2[NSpFdiff_orig>=-0.3 & NSpFdiff_orig<=-0.1])[2]

#GENFREQ
#expected
sum(GFdiff_orig>0 & GFdiff_orig<=1) #334
(PHRASE[GFdiff_orig>0 & GFdiff_orig<=1])[222]
(PHRASE_DEF[GFdiff_orig>0 & GFdiff_orig<=1])[222]
(log(ALLFREQ1)[GFdiff_orig>0 & GFdiff_orig<=1])[222]
(log(ALLFREQ2)[GFdiff_orig>0 & GFdiff_orig<=1])[222]

sum(GFdiff_orig>3 & GFdiff_orig<=4) #115
(PHRASE[GFdiff_orig>3 & GFdiff_orig<=4])[26]
(PHRASE_DEF[GFdiff_orig>3 & GFdiff_orig<=4])[26]
(log(ALLFREQ1)[GFdiff_orig>3 & GFdiff_orig<=4])[26]
(log(ALLFREQ2)[GFdiff_orig>3 & GFdiff_orig<=4])[26]

sum(GFdiff_orig>7 & GFdiff_orig<=9) #3
(PHRASE[GFdiff_orig>7 & GFdiff_orig<=9])
(PHRASE_DEF[GFdiff_orig>7 & GFdiff_orig<=9])
(log(ALLFREQ1)[GFdiff_orig>7 & GFdiff_orig<=9])
(log(ALLFREQ2)[GFdiff_orig>7 & GFdiff_orig<=9])

#unexpected
sum(GFdiff_orig>=-1 & GFdiff_orig<0) #244
(PHRASE[GFdiff_orig>=-1 & GFdiff_orig<0])[111]
(PHRASE_DEF[GFdiff_orig>=-1 & GFdiff_orig<0])[111]
(log(ALLFREQ1)[GFdiff_orig>=-1 & GFdiff_orig<0])[111]
(log(ALLFREQ2)[GFdiff_orig>=-1 & GFdiff_orig<0])[111]

sum(GFdiff_orig>=-4 & GFdiff_orig<=-3) #102
(PHRASE[GFdiff_orig>=-4 & GFdiff_orig<=-3])[22]
(PHRASE_DEF[GFdiff_orig>=-4 & GFdiff_orig<=-3])[22]
(log(ALLFREQ1)[GFdiff_orig>=-4 & GFdiff_orig<=-3])[22]
(log(ALLFREQ2)[GFdiff_orig>=-4 & GFdiff_orig<=-3])[22]

sum(GFdiff_orig>=-9 & GFdiff_orig<=-7) #3
(PHRASE[GFdiff_orig>=-9 & GFdiff_orig<=-7])
(PHRASE_DEF[GFdiff_orig>=-9 & GFdiff_orig<=-7])
(log(ALLFREQ1)[GFdiff_orig>=-9 & GFdiff_orig<=-7])
(log(ALLFREQ2)[GFdiff_orig>=-9 & GFdiff_orig<=-7])

uniquephrases<-unique(data.frame(PHRASE, PHRASE_DEF))
write.table(uniquephrases, file="../KAIST.adjadjnphrases.csv", quote=F, sep="\t")

PHRASE[ORDER=="predict"][1009:1014]
PHRASE_DEF[ORDER=="predict"][1009:1014]
RLdiff[ORDER=="predict"][1009:1014]
#1  4  2  0  2 -1

RLdiff[ORDER=="nonpredict"][1009:1014]
#-1 -4 -2  0 -2  1

#Find any examples where both orders are attested
combomods12<-paste(MOD1.ROOT, MOD2.ROOT, sep="")
combomods21<-paste(MOD2.ROOT, MOD1.ROOT, sep="")
combomods21matches<-c()
for (i in combomods12){
	if (length(which(combomods21==i))==0) {
		combomods21matches<-c(combomods21matches, 0)	
	} else {
		combomods21matches<-c(combomods21matches, which(combomods21==i)[1])	
	}
}
combomods21matches ;length(combomods21matches)
#650 adj root tokens have both orders attested

PHRASE_DEF[combomods21matches][622]
PHRASE_DEF[combomods21matches!=0][622]
PHRASE_DEF[combomods21matches][grep("red", PHRASE_DEF[combomods21matches])]

#finding a pair where "small" switches position (round small, small round)
PHRASE_DEF[combomods21matches][grep("small", PHRASE_DEF[combomods21matches], perl=T)][64]
PHRASE_DEF[combomods21matches!=0][grep("small", PHRASE_DEF[combomods21matches], perl=T)][64]
#
for (i in combomods12){
	matches<-c(matches, which(combomods21==i))
}
matches.uniq<-unique(matches)
PHRASE_DEF[matches.uniq]

combomods12.uniq<-unique(combomods12)
combomods21.uniq<-unique(combomods21)
uniq.combo.matches<-c()
for (i in combomods12.uniq){
	uniq.combo.matches <-c(uniq.combo.matches, which(combomods21.uniq==i))
}
PHRASE_DEF[uniq.combo.matches]
length(uniq.combo.matches) 
NOUN.ROOT[uniq.combo.matches]
#262/1103 uniq adj root pairs (23.7%) have both orders attested

combomods12<-paste(MOD1.ROOT, MOD2.ROOT, sep="")
combomods21<-paste(MOD2.ROOT, MOD1.ROOT, sep="")
combomods21matches<-c()
for (i in combomods12){
	if (length(which(combomods21==i))==0) {
		combomods21matches<-c(combomods21matches, 0)	
	} else {
		combomods21matches<-c(combomods21matches, which(combomods21==i)[1])	
	}
}

#PHRASE_DEF[1];PHRASE_DEF[1736]d
#Finding color ones
PHRASE_DEF[matches.uniq[grep("(white)|(red)|(blue)|(green)|(black)|(yellow)",PHRASE_DEF[matches.uniq],perl=T)]]

#Are there any where the exact same adj pair + N combo is attested in both orders? 
combomodn12<-paste(MOD1.ROOT, MOD2.ROOT, NOUN.ROOT, sep="")
combomodn21<-paste(MOD2.ROOT, MOD1.ROOT, NOUN.ROOT, sep="")
n.matches<-c()
for (i in combomodn12){
	if (length(which(combomodn21==i))==0) {
		n.matches <-c(n.matches, 0)	
	} else {
		n.matches <-c(n.matches, which(combomodn21==i)[1])	
	}
}
sum(n.matches!=0) #70
#70 adj adj n combo tokens are attested with both adj orders
sum()
PHRASE_DEF[n.matches]
PHRASE_DEF[n.matches!=0]
#how many of these are from the same file? 
sum(FILE[n.matches]==FILE[n.matches!=0]) #13 of these reverse order pairs are from the same file
#how many are the same SUBJOBJ category? 
subjobj12combo<-paste(SUBJOBJ1, SUBJOBJ2, sep="")
sum(subjobj12combo[n.matches]==subjobj12combo[n.matches!=0])  #37/70 (52.8%) are same SUBJOBJ 
#curr.kaist.file<-scan(file="./kaistcorpus_written_raw_or_literature_juvenileAndfable_CAD002.txt.new.pos2", what=character(0), sep="\n", quiet=TRUE)
line<-curr.kaist.file[grep("차갑고\t[^\t]+\t낮은", curr.kaist.file, perl=T)]
gsub("\t", " ", line, perl=T)
#[619] 매우 낮고 차가운 목소리였다  "It was a very low, cold voice"
#[620] 유의태는 더욱 차갑고 낮은 목소리로 말했다 "Yuiitae said in an even more cool low voice"