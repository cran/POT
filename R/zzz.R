".First.lib" <-
function(lib, pkg)
{
  library.dynam("POT", package = pkg, lib.loc = lib)
  return(invisible(0))
}

