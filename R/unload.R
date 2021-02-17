.onUnload <- function (libpath) {
  library.dynam.unload("match2C", libpath)
}
