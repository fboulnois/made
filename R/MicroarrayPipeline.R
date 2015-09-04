ma.pipeline <- function(config)
{
  config <- read.yaml.config(config)

  if(is.character(config$pipeline$normalization) || config$pipeline$normalization)
  {
    eset <- ma.normalize(config)
  }

  if(config$pipeline$summarization)
  {
    tf <- ma.summarize(config, eset)
  }

  if(config$pipeline$write_report)
  {
    write.report(config, eset, tf)
  }
}
