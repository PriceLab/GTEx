library(TrenaViz)
tv <- TrenaViz("TrenaProjectTemplate")
runApp(createApp(tv, port=5838))
later(function(){browseURL("http://0.0.0.0:5838")}, 2)
