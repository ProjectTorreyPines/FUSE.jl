# Stipple UI components for FUSE web GUI
# (All imports are handled in FUSEGenieFrameworkExt.jl)

# Define reactive model
@app begin
    # Case selection
    @in selected_case = "KDEMO"
    @in case_args_text = ""
    @in case_kwargs_text = ""
    @out available_cases = ["KDEMO", "ITER", "FPP"]

    # Actor selection
    @in selected_actor_edit = "ActorNoOperation"  # For parameter editing (left column)
    @in selected_actor_run = "ActorNoOperation"   # For running (right column)
    @out available_actors = ["ActorNoOperation", "ActorFluxMatcher"]

    # Parameters (dynamic, populated based on selected actor)
    @in param_values = Dict{String,String}()
    @out param_infos = Vector{Dict{String,Any}}()
    @in param_change_trigger = 0
    @out param_validation_states = Dict{String,String}()

    # Execution
    @in load_case_clicked = false
    @out is_loading_case = false
    @in run_clicked = false
    @out run_status = ""
    @out is_running = false

    # Log
    @out log_text = ""
    @in clear_log_clicked = false

    # Plot
    @out plot_image = ""

    # Session (private, not exposed to UI)
    @private session = WebGuiSession()
    @private param_values_prev = Dict{String,String}()
    @private param_values_initial = Dict{String,String}()

    # Load case when button clicked
    @onchange load_case_clicked begin
        if load_case_clicked
            is_loading_case = true
            case_sym = Symbol(selected_case)

            # Parse kwargs from text field (simple format: "key1=val1, key2=val2")
            kwargs = Dict{Symbol,Any}()
            if !isempty(case_kwargs_text)
                try
                    # Simple parsing for now - can be enhanced later
                    for pair in split(case_kwargs_text, ",")
                        if contains(pair, "=")
                            k, v = split(strip(pair), "=", limit=2)
                            kwargs[Symbol(strip(k))] = strip(v)
                        end
                    end
                catch e
                    push!(session.log_buffer, "⚠ Warning: Could not parse kwargs, using defaults")
                end
            end

            # Load the case
            success = load_case!(session, case_sym, Dict(), kwargs)

            if success && !isnothing(session.act)
                # Update parameter forms for selected actor
                actor_sym = Symbol(selected_actor)
                param_list = introspect_actor_params(session.act, actor_sym)

                # Convert ParamInfo to Dict for Stipple
                param_infos = [
                    Dict{String,Any}(
                        "name" => string(p.name),
                        "param_type" => string(p.param_type),
                        "value_type" => p.value_type,
                        "current_value" => string(p.current_value),
                        "default_value" => string(p.default_value),
                        "options" => string.(p.options),
                        "units" => p.units,
                        "description" => p.description
                    ) for p in param_list
                ]

                # Initialize param_values with current values
                param_values = Dict{String,String}(
                    string(p.name) => string(p.current_value) for p in param_list
                )
                # Initialize previous values for change detection
                param_values_prev = copy(param_values)
                # Initialize initial values for comparison
                param_values_initial = copy(param_values)
                # Initialize all validation states as pristine
                param_validation_states = Dict{String,String}(
                    string(p.name) => "pristine" for p in param_list
                )
            end

            # Update log
            log_text = join(session.log_buffer, "\n")
            is_loading_case = false
            load_case_clicked = false
        end
    end

    # Update params when actor selection changes
    @onchange selected_actor_edit begin
        if !isnothing(session.act)
            actor_sym = Symbol(selected_actor_edit)
            param_list = introspect_actor_params(session.act, actor_sym)

            # Convert ParamInfo to Dict for Stipple
            param_infos = [
                Dict{String,Any}(
                    "name" => string(p.name),
                    "param_type" => string(p.param_type),
                    "value_type" => p.value_type,
                    "current_value" => string(p.current_value),
                    "default_value" => string(p.default_value),
                    "options" => string.(p.options),
                    "units" => p.units,
                    "description" => p.description
                ) for p in param_list
            ]

            # Initialize param_values with current values
            param_values = Dict{String,String}(
                string(p.name) => string(p.current_value) for p in param_list
            )
            # Initialize previous values for change detection
            param_values_prev = copy(param_values)
            # Initialize initial values for comparison
            param_values_initial = copy(param_values)
            # Initialize all validation states as pristine
            param_validation_states = Dict{String,String}(
                string(p.name) => "pristine" for p in param_list
            )
        end
    end

    # Detect uncommitted changes as user types (turns yellow)
    @onchange param_values begin
        if !isnothing(session.act) && !isempty(param_values)
            for (param_name, param_value) in param_values
                prev_value = get(param_values_prev, param_name, "")
                initial_value = get(param_values_initial, param_name, "")
                current_state = get(param_validation_states, param_name, "pristine")

                # Only update state if it's not already invalid (preserve red on invalid)
                if current_state != "invalid"
                    if param_value != prev_value
                        # User typed something different from committed value - mark as edited (yellow)
                        param_validation_states[param_name] = "edited"
                    elseif param_value == initial_value
                        # User typed back to initial value - mark as pristine
                        param_validation_states[param_name] = "pristine"
                    end
                end
            end
        end
    end

    # Validate and commit changes when Enter is pressed
    @onchange param_change_trigger begin
        if !isnothing(session.act) && !isempty(param_values) && param_change_trigger > 0
            actor_sym = Symbol(selected_actor_edit)

            # Find which parameter(s) changed
            for (param_name, param_value) in param_values
                prev_value = get(param_values_prev, param_name, "")
                if param_value != prev_value
                    try
                        # Validate by attempting to update
                        update_param!(session.act, actor_sym, Symbol(param_name), param_value)

                        # SUCCESS: Mark as pristine (committed, no longer edited)
                        param_validation_states[param_name] = "pristine"

                        # Log the change
                        push!(session.log_buffer, "Updated act.$selected_actor_edit.$param_name = $param_value")

                    catch e
                        # FAILURE: Mark as invalid (red)
                        param_validation_states[param_name] = "invalid"
                        push!(session.log_buffer, "✗ Invalid value for $param_name: $e")
                    end
                end
            end

            # Update previous state (now committed)
            param_values_prev = copy(param_values)

            # Update the log display
            log_text = join(session.log_buffer, "\n")
        end
    end

    # Run actor when button clicked
    @onchange run_clicked begin
        if run_clicked
            is_running = true

            # Run actor (parameters are already updated in real-time)
            actor_sym = Symbol(selected_actor_run)
            success = run_actor!(session, actor_sym)

            run_status = session.last_run_status
            log_text = join(session.log_buffer, "\n")

            # Generate plot (placeholder for now)
            if success
                plot_image = plot_to_base64(session, :equilibrium)
            end

            is_running = false
            run_clicked = false
        end
    end

    # Clear log when button clicked
    @onchange clear_log_clicked begin
        if clear_log_clicked
            empty!(session.log_buffer)
            log_text = ""
            clear_log_clicked = false
        end
    end
end

# UI layout function
function ui()
    [
        Html.style("""
            .param-edited {
                background-color: #fff9c4 !important;  /* Yellow 100 */
                border-color: #fbc02d !important;       /* Yellow 700 */
            }
            .param-invalid {
                background-color: #ffcdd2 !important;  /* Red 100 */
                border-color: #c62828 !important;       /* Red 800 */
            }
        """),

        heading("FUSE Web GUI"),

        row([
            # Left panel: Controls
            cell(class="col-md-6 col-sm-12", [
                h5("Case Selection"),
                select(:selected_case; options=:available_cases, label="Select Case"),
                textfield("Kwargs (e.g., init_from=scalars)", :case_kwargs_text),
                btn("Load Case", @click(:load_case_clicked), color="primary", icon="download",
                    disable=:is_loading_case),
                Html.div(@iif(:is_loading_case), "⏳ Loading case...", class="text-grey-6 q-mt-sm"),

                separator(),

                h5("Actor Parameters"),
                select(:selected_actor_edit; options=:available_actors, label="Select Actor to Edit"),
                Html.div(class="q-pa-md", [
                    # Display parameters dynamically with inline layout
                    Html.div(
                        @recur("param in param_infos"),
                        [
                            Html.div(class="row items-center q-mb-sm", [
                                # Left column: Parameter name + help icon
                                Html.div(class="col-5", [
                                    Html.span(class="text-weight-bold", "{{ param.name }}"),
                                    Html.span(class="q-ml-xs", [
                                        icon("help_outline", size="sm", class="cursor-pointer text-grey-7",
                                            [StippleUI.tooltip(
                                                Html.div([
                                                    # Description (main content)
                                                    Html.div(
                                                        var"v-html"="param.description",
                                                        style="white-space: pre-wrap; margin-bottom: 8px;"
                                                    ),
                                                    # Metadata section
                                                    Html.div(
                                                        style="border-top: 1px solid rgba(255,255,255,0.3); padding-top: 8px; font-size: 12px; opacity: 0.9;",
                                                        [
                                                            # Type
                                                            Html.div([
                                                                Html.span("Type: ", style="font-weight: bold;"),
                                                                Html.span("{{ param.value_type }}")
                                                            ]),
                                                            # Default value
                                                            Html.div(@iif("param.default_value !== ''"), [
                                                                Html.span("Default: ", style="font-weight: bold;"),
                                                                Html.span("{{ param.default_value }}")
                                                            ]),
                                                            # Units (if present)
                                                            Html.div(@iif("param.units !== ''"), [
                                                                Html.span("Units: ", style="font-weight: bold;"),
                                                                Html.span("{{ param.units }}")
                                                            ])
                                                        ]
                                                    )
                                                ]),
                                                class="text-body2",
                                                style="font-size: 14px; max-width: 400px; padding: 12px; line-height: 1.4;"
                                            )]
                                        )
                                    ])
                                ]),

                                # Right column: Input widget
                                Html.div(class="col-7", [
                                    # Switch parameters
                                    Html.div(@iif("param.param_type === 'switch'"), [
                                        select(
                                            Symbol("param_values[param.name]"),
                                            options=Symbol("param.options"),
                                            outlined=true,
                                            dense=true,
                                            [
                                                Html.input(
                                                    type="hidden",
                                                    var"v-model"="param_values[param.name]",
                                                    var":value"="param.current_value"
                                                )
                                            ]
                                        )
                                    ]),
                                    # Entry parameters
                                    Html.div(@iif("param.param_type === 'entry'"), [
                                        Html.input(
                                            var":class"="{
                                                'q-field__native q-placeholder': true,
                                                'param-edited': param_validation_states[param.name] === 'edited',
                                                'param-invalid': param_validation_states[param.name] === 'invalid'
                                            }",
                                            style="padding: 8px 12px; border-radius: 4px; width: 100%; font-size: 14px;",
                                            var":style"="{
                                                border: param_validation_states[param.name] === 'invalid' ? '1px solid #c62828' :
                                                        param_validation_states[param.name] === 'edited' ? '1px solid #fbc02d' :
                                                        '1px solid #9e9e9e',
                                                backgroundColor: param_validation_states[param.name] === 'invalid' ? '#ffcdd2' :
                                                                param_validation_states[param.name] === 'edited' ? '#fff9c4' :
                                                                'white'
                                            }",
                                            type="text",
                                            var"v-model"="param_values[param.name]",
                                            var":placeholder"="'Default: ' + param.default_value",
                                            var"@keyup.enter"="param_change_trigger++"
                                        )
                                    ])
                                ])
                            ])
                        ]
                    ),
                    Html.p(@iif("param_infos.length === 0"), "Load a case to see parameters", class="text-grey-6")
                ])
            ]),

            # Right panel: Execution and Results
            cell(class="col-md-6 col-sm-12", [
                # Run Actor Section
                h5("Run Actor"),
                Html.div(class="q-pa-md", [
                    row([
                        cell(class="col-8", [
                            select(:selected_actor_run; options=:available_actors, label="Select Actor to Run")
                        ]),
                        cell(class="col-4", [
                            btn("Run", @click(:run_clicked), color="primary", icon="play_arrow",
                                disable=:is_running, style="width: 100%;")
                        ])
                    ]),
                    Html.div(@iif(:is_running), "⏳ Running...", class="text-grey-6 q-mt-sm")
                ]),

                separator(),

                # Status Section
                h5("Status"),
                Html.div(class="q-pa-md", [
                    Html.div(@iif("run_status === 'success'"), "✓ Actor ran successfully", class="text-positive"),
                    Html.div(@iif("run_status === 'error'"), "✗ Error running actor", class="text-negative")
                ]),

                separator(),

                # Log Section
                row([
                    cell(class="col-8", [h5("Log")]),
                    cell(class="col-4 text-right", [
                        btn("Clear Log", @click(:clear_log_clicked), color="negative", size="sm", icon="delete")
                    ])
                ]),
                Html.div(class="q-pa-md bg-grey-2", style="border: 1px solid #ccc; border-radius: 4px; max-height: 500px; overflow-y: auto; font-family: monospace; white-space: pre-wrap;", [
                    Html.pre(style="margin: 0;", "{{ log_text }}")
                ]),

                separator(),

                # Plot Section
                h5("Plot"),
                Html.div(class="q-pa-md", [
                    Html.img(src=:plot_image, @iif("plot_image !== ''"), style="max-width: 100%; height: auto;"),
                    Html.p(@iif("plot_image === ''"), "No plot available", class="text-grey-6")
                ])
            ])
        ])
    ]
end
