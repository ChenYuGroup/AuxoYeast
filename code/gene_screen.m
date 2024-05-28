
%%

model = changeRxnBounds(model,'r_1714',-1000,'l');
FBAsolution_0 = optimizeCbModel(model,'max');
miug_0=FBAsolution_0.f;
exchange=table2cell(readtable('part_temp.xlsx','VariableNamingRule','preserve'));
genes=exchange(:,1);
ID=exchange(:,4);
n = length(exchange);
result=zeros(n,1);
for i=1:n
    model_k = model;
    [model_k, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes(model_k,genes(i));
    FBAsolution_k = optimizeCbModel(model_k,'max');
    miug_k=FBAsolution_k.f;

    model_k = changeRxnBounds(model_k,ID(i),-1000,'l');%
    FBAsolution = optimizeCbModel(model_k,'max');
    miug=FBAsolution.f;

    if miug_k < 0.01*miug_0
        if miug < 0.01*miug_0
            result(i)=0;
        else
            result(i)=1;
        end
    else
        if miug < 0.01*miug_0
            result(i)=22;
        else
            result(i)=23;
        end
    end
end
