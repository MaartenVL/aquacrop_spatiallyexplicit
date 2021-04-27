function ParamStruct=readSoil4os(ParamStruct,ST,ST2,CN,texture)
ParamStruct.Soil.zSoil=ST;
ParamStruct.Soil.CN=CN;
ParamStruct.Soil.zRes=ST2;
ParamStruct.Soil.Comp.dz(:)=ST/12; % ipv 12 bv 5?
ParamStruct.Soil.Comp.dzsum=round(100*(cumsum(ParamStruct.Soil.Comp.dz)))/100; % round zorg voor te grote ST (0.7 ipv 0.697) als je diepte 3 cijfers na comma is

ParamStruct.Soil.Layer.dz=ST;
ParamStruct.Soil.Layer.Sand=texture(1);
ParamStruct.Soil.Layer.Clay=texture(2);
ParamStruct.Soil.Layer.OrgMat=texture(3);

ParamStruct;
end

