# Script for anything run on package loading

.onLoad <- function(...) {
  S7::methods_register()
}