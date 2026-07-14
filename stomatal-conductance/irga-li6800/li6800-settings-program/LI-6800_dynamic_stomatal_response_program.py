from bpdefs import EXEC, GROUP, DIALOG, IF, RETURN, ASSIGN, LOG, SETCONTROL, LOOP, Nothing, CheckBox, Text, Button, DropDown, RadioBtns, EditBox, DataDict

steps=[
# Compile and run python: EXEC(scope, [source=code] or [file=filename])
EXEC(0,file="/home/licor/resources/lib/dyn_tools.py"),
# Collection: GROUP(enabled, label)
GROUP("True","Opening Dialog",
	steps=(
		# Open a dialog: DIALOG(title=string [,sub=string] [,text=string] [,items=list_of_edit_items] [,buttons=list_of_buttons] [var=pressed_btn_name])
		DIALOG(title="'Range match & log file check'",
			sub="'Have you performed a range match? \\n Is a log file open?'",
			text="'A range match for CO2 and H2O must be run \\n every 2 hours or after a significant change in ambient temperature'",
			buttons="'Cancel','Continue'",
			var="button2"),
		IF("button2 == 'Cancel'",
			steps=(
				RETURN(),
			)
		),
		# Assign a variable to an expression: ASSIGN('varname', exp="expression" [,dlg=Nothing()])
		ASSIGN("currentplant",
			exp="'Plant ID'",
			dlg=EditBox("'Plant ID'", desc="'Enter description of plant:'")),
		# Assign a variable to an expression: ASSIGN('varname', exp="expression" [,dlg=Nothing()])
		ASSIGN("log_interval",
			exp="30",
			dlg=EditBox("'Logging interval'", units="'s'", desc="'Log every __ secs during ramp'")),
		# Collection: GROUP(enabled, label)
		GROUP("True","Last time defaults",
			steps=(
				# Compile and run python: EXEC(scope, [source=code] or [file=filename])
				EXEC(0,source="BP.setDlgFile('/home/licor/apps/dynamic/_dat_co2_continuous.json')"),
			)
		),
		# Open a dialog: DIALOG(title=string [,sub=string] [,text=string] [,items=list_of_edit_items] [,buttons=list_of_buttons] [var=pressed_btn_name])
		DIALOG(title="'Stomatal response curve'",
			items="currentplant",
			buttons="'Cancel','Continue'",
			var="button"),
		IF("button == 'Cancel'",
			steps=(
				RETURN(),
			)
		),
		# Collection: GROUP(enabled, label)
		GROUP("True","Check validity",
			steps=(
				# Compile and run python: EXEC(scope, [source=code] or [file=filename])
				EXEC(0,source="issues=checkValidity(True)"),
				IF("len(issues)",
					steps=(
						# Open a dialog: DIALOG(title=string [,sub=string] [,text=string] [,items=list_of_edit_items] [,buttons=list_of_buttons] [var=pressed_btn_name])
						DIALOG(title="'Not Ready to Run'",
							sub="'Fix the following issues, then try again:'",
							text="'\\n'.join(issues)",
							buttons="'Ok'",
							var="button"),
						RETURN(),
					)
				),
			)
		),
	)
),
# Log a remark: LOG(rem=string)
LOG(rem="currentplant"),
# Collection: GROUP(enabled, label)
GROUP("True","Collect New Data",
	steps=(
		# Compile and run python: EXEC(scope, [source=code] or [file=filename])
		EXEC(0,source="BP.setLogBuffering(15)"),
		# Collection: GROUP(enabled, label)
		GROUP("True","This sets log averaging to match log interval",
			steps=(
				# Set a control: SETCONTROL('target', 'value', 'eval' [,opt_target=''])
				SETCONTROL("LogOpts:AvgTime","True,min(15, log_interval)","int"),
			)
		),
		# Log a data record: LOG([avg='Default'] [,match='Default'] [,matchH2O='Default'] [,flr='Default'] [flash='Default'])
		LOG(match="Off",
			matchH2O="Off",
			flr="0: Nothing"),
		# Set a control: SETCONTROL('target', 'value', 'eval' [,opt_target=''])
		SETCONTROL("Qin","0","float"),
		# Loop for a duration: LOOP(dur="float" [,units='Seconds' Second|Minutes|Hours ] [,var=''] [,mininc=''])
		LOOP(dur="60",
			units="Minutes",
			mininc="30",
			steps=(
				# Log a data record: LOG([avg='Default'] [,match='Default'] [,matchH2O='Default'] [,flr='Default'] [flash='Default'])
				LOG(match="Off",
					matchH2O="Off",
					flr="0: Nothing"),
			)
		),
		# Set a control: SETCONTROL('target', 'value', 'eval' [,opt_target=''])
		SETCONTROL("Qin","1000","float"),
		# Loop for a duration: LOOP(dur="float" [,units='Seconds' Second|Minutes|Hours ] [,var=''] [,mininc=''])
		LOOP(dur="60",
			units="Minutes",
			mininc="30",
			steps=(
				# Log a data record: LOG([avg='Default'] [,match='Default'] [,matchH2O='Default'] [,flr='Default'] [flash='Default'])
				LOG(match="Off",
					matchH2O="Off",
					flr="0: Nothing"),
			)
		),
		# Set a control: SETCONTROL('target', 'value', 'eval' [,opt_target=''])
		SETCONTROL("Qin","0","float"),
		# Loop for a duration: LOOP(dur="float" [,units='Seconds' Second|Minutes|Hours ] [,var=''] [,mininc=''])
		LOOP(dur="60",
			units="Minutes",
			mininc="30",
			steps=(
				# Log a data record: LOG([avg='Default'] [,match='Default'] [,matchH2O='Default'] [,flr='Default'] [flash='Default'])
				LOG(match="Off",
					matchH2O="Off",
					flr="0: Nothing"),
			)
		),
	)
),
# Set a control: SETCONTROL('target', 'value', 'eval' [,opt_target=''])
SETCONTROL("Qin","1000","float"),
RETURN(),
]
