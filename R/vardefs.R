
##########################
# Variable Defines
# ########################
# Nunes et al, 2019
# Last updated version: 0.3.0

# Temporary directory
temp.dir = tempdir()

# Current language
#is.english = TRUE

# Verbose messages as default (FALSE in GEAP run-time)
verbose.default = TRUE

# Tags for C#<=>R connection 
tag.error = "!E!"
tag.status = "!S!"
tag.percent = "!P!"
tag.memory = "!M!"
tag.package = "!K!"

# Option tags
options(tag.error = "!E!",
        tag.status = "!S!",
        tag.percent = "!P!",
        tag.memory="!M!",
        tag.package="!K!",
        tag.output="!0!")

# Sets the API mode for GEAP
options(geap.is.api = getOption('geap.is.api', FALSE))

# "Always SVG" (no WMF) option for some of the graphics rendering
options(force.svg = TRUE)

# Local command storage (GEAP-only use)
.command.buffer = raw(4096L)
.command.last.return = NULL
.last.serialized = NULL

