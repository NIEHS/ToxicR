# Copyright 2021  NIEHS <matt.wheeler@nih.gov>
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software
# and associated documentation files (the Software), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge, publish, distribute,
# sublicense, and/or sell copies of the Software, and to permit persons to whom the Software
# is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies
# or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED AS IS, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
# CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
# OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
.onAttach <- function(libname, pkgname) {
  msg <- "

  _______               _____
 |__   __|      \U1F913     |  __ \\
    | | _____  ___  ___| |__) |
    | |/ _ \\ \\/ / |/ __|  _  /
    | | (_) >  <| | (__| | \\ \\
    |_|\\___/_/\\_\\_|\\___|_|  \\_\\
                           23.10.1.1.1
    ___
    | |
    / \\                   ____()()
   /\U2620\U2620\U2620\\                 /       xx
  (_____)          `~~~~~\\_;m__m.__>o


THE SOFTWARE IS PROVIDED AS IS, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"

  # notify the option to update? (5 possible cases)
  # 1. error in getting the current package version -> yes
  # 2. error in getting the github package version -> yes
  # 3. current package version < github package version -> yes
  # 4. current package version > github package version -> no
  # 5. current package version == github package version -> no
  # in short, notify the option to update unless the version numbers match
  # get version of the currently installed package
  current_pkg_version <- tryCatch(
    as.character(utils::packageVersion("ToxicR")),
    error = function(e) "unknown")
  # github url
  github_url <- paste0(
    "https://raw.githubusercontent.com/SciomeLLC/ToxicR/",
    "main/DESCRIPTION")
  # get github description or handle errors
  github_pkg_desc <- tryCatch(
    readLines(github_url),
    warning = function(w) {"github_desc_read_fail"},
    error = function(e) {"github_desc_read_fail"})
  # get the version number of github version
  if (identical(github_pkg_desc, "github_desc_read_fail")) {
    github_pkg_version <- "unknown"
  } else {
    github_pkg_version <- gsub(
      ".*ersion: ", "", github_pkg_desc[
        grep("ersion", github_pkg_desc)][1])
  }
  # compare versions
  compare_version_result <- tryCatch(
    utils::compareVersion(
      current_pkg_version, github_pkg_version),
    warning = function(w) {999}, # 999 indicates no need for update
    error = function(e) {999})
  # skip update for case 5
  if (current_pkg_version != "unknown" &
      github_pkg_version != "unknown" &
      compare_version_result == 0) {
    startup_message <- paste0(
      "Package attached: ToxicR v", current_pkg_version,
      " (same as the most recent version available through GitHub).")
  } else if (
    # skip update for case 4
    current_pkg_version != "unknown" &
    github_pkg_version != "unknown" &
    compare_version_result > 0) {
    startup_message <- paste0(
      "Package attached: ToxicR v", current_pkg_version,
      " (probably the most recent version available through GitHub).")
  } else {
    # simply notify of the OPTION to update the package
    # this is simply a notification of the option to update,
    # rather than a recommendation to update
    startup_message <- paste0(
      "Package attached: ToxicR v", current_pkg_version,
      "; Most recent version available on GitHub: v", github_pkg_version,
      "\n\nYou have an option to update the package ",
      "with the function `update_ToxicR()`. ",
      "If you do so, make sure to restart R.\n\n")
  }
  packageStartupMessage(msg)
  packageStartupMessage(startup_message)
  .set_threads(2);
  open_mp_message <- "OpenMP threads set to 2.\n"
  packageStartupMessage(open_mp_message)
}
