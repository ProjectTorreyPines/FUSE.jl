# Efficient ClaudeCode development in Julia with RemoteREPL.jl and Revise.jl

## Purpose
Use RemoteREPL.jl and Revise.jl for debugging Julia packages with long startup times by maintaining a permanently running Julia session. This avoids repeated package loading and compilation overhead.

## Installation
```julia
using Pkg
Pkg.add("RemoteREPL")
Pkg.add("Revise")
```

## Basic RemoteREPL Setup

### Start RemoteREPL Server
In your main Julia session (this will be the persistent session):

```julia
using RemoteREPL, Sockets
@async serve_repl(Sockets.localhost, 9999)
```

### Connect to Remote Session and Command Execution
From another Julia process or session:

```julia
using RemoteREPL, Sockets
con = connect_remote(Sockets.localhost, 9999)

# Simple commands
result = remotecmd(con, "2 + 2")  # Returns: 4

# Variable creation and manipulation
remotecmd(con, "x = 42")
result = remotecmd(con, "x * 2")  # Returns: 84

# Multi-line commands
remotecmd(con, """
for i in 1:3
    println("Line ", i)
end
"finished"
""")
```

### Output and Printing Behavior
**Important**: RemoteREPL handles output differently:

- **Return values**: Only the last expression is returned by `remotecmd()`
- **println() output**: Goes to server console, NOT captured in return value
- **Logging macros**: `@info`, `@warn`, `@error` ARE captured and included in result
- **display() output**: IS captured and included in result

```julia
# This println goes to server console, not returned
remotecmd(con, "println(\"Hello\"); 123")  # Returns: 123

# Logging is captured
remotecmd(con, "@info \"This appears\"; 456")  # Returns: "[ Info: This appears\n456"

# Display output is captured  
remotecmd(con, "display([1,2,3]); 789")  # Returns: "3-element Vector{Int64}:\n 1\n 2\n 3\n789"
```

### Error Handling

RemoteREPL returns errors as strings rather than throwing them locally:

```julia
# Syntax errors
result = remotecmd(con, "x = 1 +")  # Returns detailed ParseError string

# Runtime errors  
result = remotecmd(con, "undefined_variable")  # Returns UndefVarError string

# Method errors
result = remotecmd(con, "sqrt(\"hello\")")  # Returns MethodError with suggestions
```

The server continues running even after errors occur.

## Iterative development in ClaudeCode using RemoteREPL.jl and Revise.jl

This approach can reduce ClaudeCode development cycle time from minutes to seconds for packages with heavy loading/precompiling dependencies such as FUSE.

1. Let ClaudeCode start a separate process with a permanent session on a specific port
2. Let ClaudeCode make code changes to package files
3. Have ClaudeCode `connect_remote` and then run tests remotely: `remotecmd(con, "include(\"test/runtests.jl\")")`
4. Debug: `remotecmd(con, "@show problematic_variable")`
5. Repeat from point 2

NOTE:
1. **Start server early**: Begin your RemoteREPL server at the start of your work session
2. **Use Revise**: Always load Revise.jl in the server for automatic code reloading
3. **Test incrementally**: Run small tests frequently rather than full test suites
4. **Monitor output**: Remember that `println()` goes to server console, use `@info` for captured output
5. **Handle errors gracefully**: Check return strings for error messages

### 1. Start Persistent Session
```julia
# Terminal 1: Start server with your slow packages
using RemoteREPL, Revise
@async serve_repl(Sockets.localhost, 9999)
```

### 2. Connect from Development Environment
```julia
# Terminal 2: Connect and work
using RemoteREPL, Sockets
con = connect_remote(Sockets.localhost, 9999)

# Run commands in persistent session
remotecmd(con, "using FUSE\")")
remotecmd(con, "ini, act = FUSE.case_parameters(:ITER; init_from=:ods)")
remotecmd(con, "dd = FUSE.init(ini, act)")
```

### 3. Connect from Development Environment
Variables and state persist between connections, even after disconnect/reconnect
```julia
# Terminal 2: Re-Connect and continue to work
using RemoteREPL, Sockets
con = connect_remote(Sockets.localhost, 9999)

remotecmd(con, "dd")  # Data still available
```
