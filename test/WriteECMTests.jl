using PhysiCellECMCreator

createICECMXMLTemplate("test")

config_dict = Dict("x_min" => -400.0, "x_max" => 400.0, "y_min" => -400.0, "y_max" => 400.0, "dx" => 40.0, "dy" => 40.0)
generateICECM("test/ecm.xml", "test.csv", config_dict)