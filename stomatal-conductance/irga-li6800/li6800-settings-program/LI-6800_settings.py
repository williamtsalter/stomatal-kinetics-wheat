from bpdefs import SETCONTROL, Nothing, CheckBox, Text, Button, DropDown, RadioBtns, EditBox, DataDict

steps=[
# Set a control: SETCONTROL('target', 'value', 'eval' [,opt_target=''])
SETCONTROL("Fan_rpm","10000","float"),
# Set a control: SETCONTROL('target', 'value', 'eval' [,opt_target=''])
SETCONTROL("Flow","600","float"),
# Set a control: SETCONTROL('target', 'value', 'eval' [,opt_target=''])
SETCONTROL("Pressure","0.05","float"),
# Set a control: SETCONTROL('target', 'value', 'eval' [,opt_target=''])
SETCONTROL("Color_All","'r90B40'","string"),
# Set a control: SETCONTROL('target', 'value', 'eval' [,opt_target=''])
SETCONTROL("Q_All","1000","float"),
# Set a control: SETCONTROL('target', 'value', 'eval' [,opt_target=''])
SETCONTROL("Txchg","22","float"),
# Set a control: SETCONTROL('target', 'value', 'eval' [,opt_target=''])
SETCONTROL("FlrOpt:Action","0: Nothing",""),
# Set a control: SETCONTROL('target', 'value', 'eval' [,opt_target=''])
SETCONTROL("Meas:Recording","False",""),
# Set a control: SETCONTROL('target', 'value', 'eval' [,opt_target=''])
SETCONTROL("Meas:Modulation","Off",""),
# Set a control: SETCONTROL('target', 'value', 'eval' [,opt_target=''])
SETCONTROL("CO2_r","400","float"),
# Set a control: SETCONTROL('target', 'value', 'eval' [,opt_target=''])
SETCONTROL("VPD_leaf","1.6","float"),
]
