model = readCbModel('yeast9.0.xml');
model = buildRxnGeneMat(model); 
model = changeRxnsModel(model);
model = changeRxnBounds(model,'r_1714',-1000,'l');
FBAsolution = optimizeCbModel(model,'max');
miug=FBAsolution.f;
genes=model.genes;
geneNames=model.geneNames;
%
exRxns = {};
exRxnNames = {};
exMetNames = {};
for i = 1:size(model.S, 2)
    if length(find(model.S(:, i))) == 1
        exRxns = [exRxns;model.rxns(i)];
        exMetNames = [exMetNames; model.metNames(find(model.S(:, i)))];
        exRxnNames = [exRxnNames;model.rxnNames(i)];
    end
end
%
result = {};
for i=1:length(genes)
    model_k = model;
    [model_k, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes(model_k,genes(i));
    FBAsolution_k = optimizeCbModel(model_k,'max');
    miug_k=FBAsolution_k.f;
    for j=1:length(exRxns)
        model_a = model_k;
        model_a = changeRxnBounds(model_a,exRxns(j),-1000,'l');
        FBAsolution_a = optimizeCbModel(model_a,'max');
        miug_a=FBAsolution_a.f;
        if miug_k < 0.01*miug && miug_a > 0.01*miug
            tmp = [geneNames(i),exMetNames(j)];
            result = [result;tmp];
        end
    end
end
writecell(result, 'g_sall.xlsx');