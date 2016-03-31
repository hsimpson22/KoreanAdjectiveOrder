KoreanAdjectiveOrder

This repo contains the final versions of my article analyzing adjective order in Korean (currently under review). The dataset used for this study is KAIST Korean corpus (http://semanticweb.kaist.ac.kr/home/index.php/KAIST_Corpus). The included morphological analyzer, trained on the KAIST corpus, was used to parse the full raw corpus (approximately 30 million words).  

Publishable.R contains the accompanying R code I wrote to process this data a few years ago. WARNING: it is VERY UGLY, but it works :)

Publishable.R does the following : 
- reads in a .csv file containing the matches to a regex I wrote to extract adjective-adjective-noun sequences in the morphologically-tagged KAIST corpus. 
- reads in another .csv file containing annotations from native speaker research assistant for those adj-adj-n sequences 
- matches the annotations up with original set of sequences
- fixes MANY one-off bad data instances (mostly due to unexpected punctuation marks that caused an error in the parser and/or regex)
- computes generalized linear regression model predicting the order of the two adjectives based on the difference between their values for 7 numeric/ordinal variables (root length, semantic openness, independence from comparison, subjectivity of adjective quality, affective load, noun-specific frequency) and one categorical variable (construction type). The details of these variables are explained in the accompanying paper. 
- evaluates the classification accuracy of the model, visualizes the model results (with library effects), and performs some post-hoc analysis
