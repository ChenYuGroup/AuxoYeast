  function model = auxoCurate(model)

% AAT2-aspartate
model = changeRxnBounds(model,'r_0217',0,'b');

% ALD2 and ALD3-pantothenic acid
g1=find(ismember(model.genes,'YMR170C'));
g2=find(ismember(model.genes,'YMR169C'));
model.rules(find(ismember(model.rxns,'r_0172'))) = {['x(' num2str(g1) ') | x(' num2str(g2) ')']};

% CHO2-choline
g1=find(ismember(model.genes,'YGR157W'));
model.rules(find(ismember(model.rxns,'r_2488'))) = {['x(' num2str(g1) ')']};
model.rules(find(ismember(model.rxns,'r_2489'))) = {['x(' num2str(g1) ')']};
model.rules(find(ismember(model.rxns,'r_2490'))) = {['x(' num2str(g1) ')']};
model.rules(find(ismember(model.rxns,'r_2491'))) = {['x(' num2str(g1) ')']};
model.rules(find(ismember(model.rxns,'r_2492'))) = {['x(' num2str(g1) ')']};
model.rules(find(ismember(model.rxns,'r_2493'))) = {['x(' num2str(g1) ')']};
model.rules(find(ismember(model.rxns,'r_2494'))) = {['x(' num2str(g1) ')']};
model.rules(find(ismember(model.rxns,'r_2495'))) = {['x(' num2str(g1) ')']};

% CYS3-cysteine & CYS4-cysteine
model = changeRxnBounds(model,'r_4703',0,'b');
model = changeRxnBounds(model,'r_0312',0,'b');

% heme a exchange ERG13-ergosterol & HEM1-heme...
[model, ~] = addExchangeRxn(model,'s_3714[c]', 0, 1000);%add heme a exchange
model.rxns(end) = {'r_temp1'};
model.rxnNames(end) = {'heme a exchange'};
% ERG10-ergosterol
model = changeRxnBounds(model,'r_0559',0,'b'); %add heme a

% GSH1-glutathione
S_current = model.S;
r = find(ismember(model.rxns,'r_4598'));%cofactor
j = find(ismember(model.mets,'s_0750[c]'));%glutathione
S_current(j, r) = -1e-06;

% GFA1-D-Glucosamine
r = find(ismember(model.rxns,'r_4048'));%
j1 = find(ismember(model.mets,'s_0509[c]'));%chitin
j2 = find(ismember(model.mets,'s_0001[ce]'));%(1->3)-beta-D-glucan
j3 = find(ismember(model.mets,'s_0004[ce]'));%(1->6)-beta-D-glucan
j4 = find(ismember(model.mets,'s_1107[c]'));%mannan
j5 = find(ismember(model.mets,'s_0773[c]'));%glycogen
j6 = find(ismember(model.mets,'s_1520[c]'));%trehalose
S_current(j1, r) = -0.02361;
S_current(j2, r) = -0.73914;
S_current(j3, r) = -0.24696;
S_current(j4, r) = -0.70204;
S_current(j5, r) = -0.35689;
S_current(j6, r) = -0.13655;
model.S = S_current;

g1=find(ismember(model.genes,'YKL104C'));
model.rules(find(ismember(model.rxns,'r_0477'))) = {['x(' num2str(g1) ')']};

% MET3-methionine
model = changeRxnBounds(model,'r_1026',0,'b');

% MET13-methionine
g1=find(ismember(model.genes,'YGL125W'));
g2=find(ismember(model.genes,'YPL023C'));
model.rules(find(ismember(model.rxns,'r_0080'))) = {['(x(' num2str(g1) ') & x(' num2str(g2) ')) | x(' num2str(g1) ')']};

% MET17-methionine
model = changeRxnBounds(model,'r_0815',0,'l');

% THI4-thiamine & HOM2-methionine and threonine...
model = changeRxnBounds(model,'r_2070',0,'b');
model = changeRxnBounds(model,'r_2071',0,'b');
model = addReaction(model,'r_temp2','reactionName','HET-P synthase (thiazole synthase)', ...
                     'reactionFormula','s_3910[c] + s_0803[c] -> s_0423[c] + s_0293[c] + s_0456[c] + s_0794[c]', ...
                     'geneRule', 'YGR144W');
model = addReaction(model,'r_temp3','reactionName','adenylated thiazole synthase', ...
                     'reactionFormula','s_1003[c] + s_1198[c] + s_0841[c] -> s_1216[c] + s_3910[c] + 3 s_0803[c] + s_0794[c]', ...
                     'geneRule', 'YGR144W');

% URA2-uracil
g1=find(ismember(model.genes,'YOR303W'));
g2=find(ismember(model.genes,'YJR109C'));
g3=find(ismember(model.genes,'YJL130C'));
model.rules(find(ismember(model.rxns,'r_0250'))) = {['(x(' num2str(g1) ') & x(' num2str(g2) ')) | x(' num2str(g3) ')']};
