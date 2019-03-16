cat('Loading proto-package ....\n')
script_dir <- dirname(sys.frame(1)[['ofile']])
tool_scripts <- list.files(path = script_dir, pattern = '.R')
tool_scripts <- file.path(script_dir, tool_scripts[tool_scripts != 'load.R'])
for (tool_script in tool_scripts) {
  print(tool_script)
  source(tool_script)
}