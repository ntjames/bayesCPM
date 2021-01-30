## made with bookdown: https://bookdown.org/yihui/bookdown/

# define output directory
sinfo <- Sys.info()
if(sinfo['sysname']=='Windows' & sinfo['user']=='nj115'){
  bayes_cpm_pap_dir <- file.path("C:/Users/nj115/Dropbox/njames/school/PhD/orise_ra/bayes_cpm/paper_arxiv")
} else if (sinfo['sysname']=='Linux' & sinfo['user']=='nathan') {
  bayes_cpm_pap_dir <- file.path("/home/nathan/Dropbox/njames/school/PhD/orise_ra/bayes_cpm/paper_arxiv")
}
  
out_dir <- file.path(bayes_cpm_pap_dir,"output")

# define working directory
setwd(bayes_cpm_pap_dir)

## html gitbook
#bookdown::render_book("index.Rmd", "bookdown::gitbook", output_dir = out_dir)

## pdf
bookdown::render_book(file.path(bayes_cpm_pap_dir,"index.Rmd"), "bookdown::pdf_book", output_dir = out_dir)

## word doc
#bookdown::render_book("index.Rmd", "bookdown::word_document2", output_dir = out_dir)

#for citations
#https://cran.r-project.org/web/packages/citr/readme/README.html

