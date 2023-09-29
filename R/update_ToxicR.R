#' Update the package 'ToxicR'
#'
#' @title update_ToxicR - downloads and tries to install the latest version of the package from github
#' Updates the current package 'ToxicR' by installing the
#' most recent version of the package from GitHub
#' This function requires installing Package 'remotes' v2.4.2
#' (or possibly a higher version) by Csardi et al. (2021),
#' <https://cran.r-project.org/package=remotes>
#'
#' @param force logical. If \code{force = TRUE}, force installing the
#' update. If \code{force = FALSE}, do not force installing the update.
#' By default, \code{force = TRUE}.
#' @param upgrade_other_pkg input for the \code{upgrade} argument to
#' be passed on to \code{remotes::install_github}.
#' One of "default", "ask", "always", "never", TRUE, or FALSE.
#' "default" respects the value of the R_REMOTES_UPGRADE environment
#' variable if set, and falls back to "ask" if unset.
#' "ask" prompts the user for which out of date packages to upgrade.
#' For non-interactive sessions "ask" is equivalent to "always".
#' TRUE and FALSE correspond to "always" and "never" respectively.
#' By default, \code{upgrade_other_pkg = FALSE}.
#' @param confirm logical. If \code{confirm = TRUE}, the user will
#' need to confirm the update. If \code{confirm = FALSE}, the confirmation
#' step will be skipped. By default, \code{confirm = TRUE}.
#' @return there will be no output from this function. Rather, executing
#' this function will update the current 'ToxicR' package by installing
#' the most recent version of the package from GitHub.
#' @examples
#' \dontrun{
#' if (interactive()) {update_ToxicR()}
#' }
#'
#' @export
update_ToxicR <- function(force = TRUE, upgrade_other_pkg = FALSE, confirm = TRUE) {
    # 6 possible cases
    # 1. error in getting the current package version -> yes
    # 2. error in getting the github package version -> yes
    # 3. current package version < github package version -> yes
    # 4. current package version > github package version -> no
    # 5. current package version == github package version -> no
    # 6. the user forces the update -> yes
    # in short, notify the option to update unless the version numbers match
    if (force == TRUE) {
        # unload the package ToxicR
        while ("package:ToxicR" %in% search()) {
            unloadNamespace("ToxicR")
        }
        remotes::install_github(
        "SciomeLLC/ToxicR", force = force, upgrade = upgrade_other_pkg)
        # attach the package
        ToxicR::prep("ToxicR", silent_if_successful = TRUE)
    } else {
        # get version of the currently installed package
    current_pkg_version <- tryCatch(
      as.character(utils::packageVersion("kim")),
      error = function(e) "unknown")
    # github url
    github_url <- paste0(
      "https://raw.githubusercontent.com/SciomeLLC/",
      "ToxicR/github-actions-build/DESCRIPTION")
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
          grep("ersion", github_pkg_desc)])
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
      message(paste0(
        "Current version of 'ToxicR': v", current_pkg_version,
        " (same as the most recent version available through GitHub)."))
    } else if (
    # skip update for case 4
    current_pkg_version != "unknown" &
    github_pkg_version != "unknown" &
    compare_version_result > 0) {
    message(paste0(
        "Current version of 'ToxicR': v", current_pkg_version,
        " (probably the most recent version available through GitHub)."))
    } else {
        # confirm update
        if (confirm == TRUE) {
            # ask the user to confirm
            user_reply <- utils::menu(
            c("Yes.", "No."),
            title = "\nDo you want to try to update the package 'ToxicR'?")
        } else {
            # if not asked, assume the user wants to update
            user_reply <- 1
        }
        # update if user wants the update
        if (user_reply == 1) {
            # unload the package kim
            while ("package:ToxicR" %in% search()) {
                unloadNamespace("ToxicR")
            }
            remotes::install_github(
            "SciomeLLC/ToxicR", force = force, upgrade = upgrade_other_pkg)
            # attach the package
            prep("ToxicR", silent_if_successful = TRUE)
            }
        }
    }
}

#' Prepare package(s) for use
#'
#' Installs, loads, and attaches package(s). If package(s) are not
#' installed, installs them prior to loading and attaching.
#'
#' @title prep
#' @param ... names of packages to load and attach, separated by commas,
#' e.g., \code{"ggplot2", data.table}. The input can be any number
#' of packages, whose names may or may not be wrapped in quotes.
#' @param pkg_names_as_object logical. If \code{pkg_names_as_object = TRUE},
#' the input will be evaluated as one object containing package names.
#' If \code{pkg_names_as_object = FALSE}, the input will be considered
#' as literal packages names (default = FALSE).
#' @param silent_if_successful logical. If \code{silent_if_successful = TRUE},
#' no message will be printed if preparation of package(s) is successful.
#' If \code{silent_if_successful = FALSE}, a message indicating which
#' package(s) were successfully loaded and attached will be printed
#' (default = FALSE).
#' @param silent_load_pkgs a character vector indicating names of
#' packages to load silently (i.e., suppress messages that get printed
#' when loading the packaged). By default, \code{silent_load_pkgs = NULL}
#' @return there will be no output from this function. Rather, packages
#' given as inputs to the function will be installed, loaded, and attached.
#' @examples
#' \donttest{
#' prep(data.table)
#' prep("data.table", silent_if_successful = TRUE)
#' prep("base", utils, ggplot2, "data.table")
#' pkgs <- c("ggplot2", "data.table")
#' prep(pkgs, pkg_names_as_object = TRUE)
#' prep("data.table", silent_load_pkgs = "data.table")
#' }
#'
prep <- function(..., pkg_names_as_object = FALSE, silent_if_successful = FALSE, silent_load_pkgs = FALSE) {
    arg_list <- as.list(match.call(expand.dots = FALSE))[["..."]]
    if ("..." %in% names(arg_list)) {
        arg_list <- arg_list[["..."]]
    }
    if (typeof(arg_list) == "language") {
        arg_list <- eval(arg_list)
    }

    # names of packages
    # check whether an object containing package names was entered
    # (e.g., a vector of package names)
    if (pkg_names_as_object == TRUE) {
        if (length(arg_list) != 1) {
            stop(paste0(
                "If pkg_names_as_object = TRUE, exactly one object containing ",
                "package names can be entered.\n",
                "You have entered ", length(arg_list), " object(s)."))
        } else {
            pkg_names <- eval.parent(arg_list[[1]])
        }
    } else {
        pkg_names <- lapply(arg_list, as.character)
    }
    if (length(pkg_names) < 1) {
        stop("Please enter name(s) of package(s) to prepare for use")
    }
    # unlist package names
    entered_pkg <- unlist(pkg_names)
    # done dealing with inputs

    # list of currently installed packages
    installed_packages_1 <- rownames(utils::installed.packages())

    pkg_to_install <- entered_pkg[
        which(!entered_pkg %in% installed_packages_1)
    ]
    if (length(pkg_to_install) > 0) {
        # install packages
        utils::install.packages(pkg_to_install)
    }
    # list of installed packages after installing entered packages
    installed_packages_2 <- rownames(utils::installed.packages())
    # newly installed packages
    newly_installed_pkg <- setdiff(
        installed_packages_2, installed_packages_1
    )
    # packages that failed to install
    pkg_that_failed_to_install <- setdiff(
        pkg_to_install, installed_packages_2
    )
    # check which packages can be loaded (and attached)
    pkg_to_check_if_loaded <- entered_pkg[
        which(entered_pkg %in% installed_packages_2)
    ]
    if (length(pkg_to_check_if_loaded) > 0) {
        # currently loaded packages
        pkg_loaded <- .packages()
        pkg_to_load <- pkg_to_check_if_loaded[
            which(!pkg_to_check_if_loaded %in% pkg_loaded)
        ]
        if (length(pkg_to_load) > 0) {
            # silently load packages
            if (length(silent_load_pkgs) > 0) {
                # split pkg_to_load into two groups: (1) packages to be loaded
                # silently and (2) packages to be loaded not silently
                pkg_to_actually_load_silently <- intersect(
                    silent_load_pkgs, pkg_to_load)
                    pkg_to_load_not_silently <- setdiff(
                    pkg_to_load, silent_load_pkgs)
                # is any of the packages to be loaded silently in pkg_to_load?
                if (length(pkg_to_actually_load_silently) > 0) {
                    invisible(lapply(pkg_to_actually_load_silently, function(x) {
                        suppressMessages(library(x, character.only = TRUE))
                    }))
                if (length(pkg_to_load_not_silently) > 0) {
                    # load and attach packages
                    invisible(lapply(
                    pkg_to_load_not_silently, library, character.only = TRUE))
                }
            } else {
                stop(paste0(
                    "None of the packages given for the argument silent_load_pkgs ",
                    "is in the packages given for the argument pkg_to_load."))
            }
        } else {
            # load and attach packages
            invisible(lapply(pkg_to_load, library, character.only = TRUE))
        }
    }
  }
    # packages that failed to install
    if (length(pkg_that_failed_to_install) > 0) {
        message(paste0(
        "The following package(s) failed to install: ",
        paste(pkg_that_failed_to_install, collapse = ", ")))
    }
    # packages that failed to load
    pkg_that_failed_to_load <- setdiff(
        entered_pkg, rownames(utils::installed.packages())
    )
    if (length(pkg_that_failed_to_load) > 0) {
        message(paste0(
            "The following package(s) failed to load: ",
            paste(pkg_that_failed_to_load, collapse = ", ")))
    }
    # newly installed packages
    if (length(newly_installed_pkg) > 0) {
        message(paste0(
            "The following package(s) were newly installed: ",
            paste(newly_installed_pkg, collapse = ", ")))
    }
    # packages that were loaded
    loaded_pkg_final <- intersect(
        entered_pkg, rownames(utils::installed.packages())
    )
    if (length(loaded_pkg_final) > 0) {
        if (silent_if_successful == FALSE) {
            message(paste0(
                "The following package(s) were successfully loaded: ",
                paste(loaded_pkg_final, collapse = ", ")))
        }
    }
}