using PhysiCellECMCreator

config_dict = Dict("x_min" => -400.0, "x_max" => 400.0, "y_min" => -400.0, "y_max" => 400.0, "dx" => 40.0, "dy" => 40.0)

createICECMXMLTemplate("test")
generateICECM("test/ecm.xml", "test.csv", config_dict)

createICECMXMLTemplate("test_monolayer"; monolayer=true)
generateICECM("test_monolayer/ecm.xml", "test_monolayer.csv", config_dict)
