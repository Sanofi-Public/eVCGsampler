# Run this manually after editing www/sass/*.scss - the app loads the
# compiled www/css/app.min.css, not the SCSS sources.
#
#   source("R/sass.R")
#
# Guarded so that sourcing R/ at app start does not recompile the CSS on
# every launch (and does not fail where sass is not installed).
if (interactive() && identical(Sys.getenv("VCG_BUILD_CSS"), "1")) {
    sass::sass(
        input = list(
            sass::sass_file("www/sass/app.scss")
        ),
        output = "www/css/app.min.css",
        options = sass::sass_options(output_style = "compressed")
    )
}
