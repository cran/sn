.onAttach <- function(library, pkg)
{
  # require("stats4")  
  # require("methods")
  # require("mnormt")
  # require("numDeriv")
  if(interactive())
  {
    meta <- packageDescription("sn")
    packageStartupMessage(
         "Package 'sn', ", meta$Version, " (", meta$Date, "). ",
         "Type 'help(SN)' for summary information.")
  }
  invisible()
}