require(shiny)
folder_address = 'C:\\efs_copy\\shiny\\ccleExplorerFull'

x <- system("ipconfig", intern=TRUE)
z <- x[grep("IPv4", x)]
ip <- gsub(".*? ([[:digit:]])", "\\1", z)
print(paste0("the Shiny Web application runs on: http://", ip, ":1234/"))

runApp(folder_address, launch.browser=FALSE, port = 1234, host = ip)
