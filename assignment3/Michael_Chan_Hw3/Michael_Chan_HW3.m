function[result,result15,resultF,resultF15] = Michael_Chan_HW3

%% read and collect the data
[im, person, number, subset] = readFaceImages('faces');

% take each 50x50 and reshape into 2500x1 vector
for i = 1:numel(im)
    im{i} = reshape(im{i},[2500,1]);
end

%% group images based on subsets
set1 = find(subset==1);
subset1 = [];
for element = 1:numel(set1)
    subset1 = horzcat(subset1,im{set1(element)});
    personlist1(element) = person(set1(element));
end

set2 = find(subset==2);
subset2 = [];
for element = 1:numel(set2)
    subset2 = horzcat(subset2,im{set2(element)});
    personlist2(element) = person(set2(element));
end
set3 = find(subset==3);
subset3 = [];
for element = 1:numel(set3)
    subset3 = horzcat(subset3,im{set3(element)});
    personlist3(element) = person(set3(element));
end
set4 = find(subset==4);
subset4 = [];
for element = 1:numel(set4)
    subset4 = horzcat(subset4,im{set4(element)});
    personlist4(element) = person(set4(element));
end
set5 = find(subset==5);
subset5 = [];
for element = 1:numel(set5)
    subset5 = horzcat(subset5,im{set5(element)});
    personlist5(element) = person(set5(element));
end

%% put together subset 1 and 5 for later training set
subset15 = horzcat(subset1,subset5);
personlist15 = horzcat(personlist1,personlist5);



%% eigenfaces

%% eigenfaces trained on subset 1
[ureduce,utot] = eigenfaces(subset1);
trainSet = applyUreduce(ureduce,subset1,utot);

testSet1 = applyUreduce(ureduce,subset1,utot);
testSet2 = applyUreduce(ureduce,subset2,utot);
testSet3 = applyUreduce(ureduce,subset3,utot);
testSet4 = applyUreduce(ureduce,subset4,utot);
testSet5 = applyUreduce(ureduce,subset5,utot);


% %%recontruction plotting
% for k=1:1:5
%     w = testSet1{k};
%     u = ureduce(:,k);
% reconstruct{k} = utot+w(k)*u;
% subplot(2,5,k)
% colormap gray
% imshow(reshape(subset1(:,k),50,50));
% subplot(2,5,k+5)
% imagesc(reshape(reconstruct{k},50,50))
% end

result(1) = nearestNeighbors(trainSet,testSet1,personlist1,personlist1);
result(2) = nearestNeighbors(trainSet,testSet2,personlist1,personlist2);
result(3) = nearestNeighbors(trainSet,testSet3,personlist1,personlist3);
result(4) = nearestNeighbors(trainSet,testSet4,personlist1,personlist4);
result(5) = nearestNeighbors(trainSet,testSet5,personlist1,personlist5);

%% eigenfaces trained on subset 1,5

[ureduce15, utot15] = eigenfaces(subset15);
trainSet15 = applyUreduce(ureduce15,subset15,utot15);

testSet115 = applyUreduce(ureduce15,subset1,utot15);
testSet215 = applyUreduce(ureduce15,subset2,utot15);
testSet315 = applyUreduce(ureduce15,subset3,utot15);
testSet415 = applyUreduce(ureduce15,subset4,utot15);
testSet515 = applyUreduce(ureduce15,subset5,utot15);

result15(1) = nearestNeighbors(trainSet15,testSet115,personlist15,personlist1);
result15(2) = nearestNeighbors(trainSet15,testSet215,personlist15,personlist2);
result15(3) = nearestNeighbors(trainSet15,testSet315,personlist15,personlist3);
result15(4) = nearestNeighbors(trainSet15,testSet415,personlist15,personlist4);
result15(5) = nearestNeighbors(trainSet15,testSet515,personlist15,personlist5);

%% fisherfaces

%% fisherfaces trained on subset1

[ureduceFt,u]= fisherFaces(subset1,personlist1);

trainSetFish = applyUreduce(ureduceFt,subset1,u);
testSet1Fish = applyUreduce(ureduceFt,subset1,u);
testSet2Fish = applyUreduce(ureduceFt,subset2,u);
testSet3Fish = applyUreduce(ureduceFt,subset3,u);
testSet4Fish = applyUreduce(ureduceFt,subset4,u);
testSet5Fish = applyUreduce(ureduceFt,subset5,u);

resultF(1) = nearestNeighbors(trainSetFish,testSet1Fish,personlist1,personlist1);
resultF(2) = nearestNeighbors(trainSetFish,testSet2Fish,personlist1,personlist2);
resultF(3) = nearestNeighbors(trainSetFish,testSet3Fish,personlist1,personlist3);
resultF(4) = nearestNeighbors(trainSetFish,testSet4Fish,personlist1,personlist4);
resultF(5) = nearestNeighbors(trainSetFish,testSet5Fish,personlist1,personlist5);

%% fisherfaces trained on subset 1,5

[ureduceF15,u15]= fisherFaces(subset15,personlist15);
trainSetFish15 = applyUreduce(ureduceF15,subset15,u15);

testSet1Fish15 = applyUreduce(ureduceF15,subset1,u15);
testSet2Fish15 = applyUreduce(ureduceF15,subset2,u15);
testSet3Fish15 = applyUreduce(ureduceF15,subset3,u15);
testSet4Fish15 = applyUreduce(ureduceF15,subset4,u15);
testSet5Fish15 = applyUreduce(ureduceF15,subset5,u15);

resultF15(1) = nearestNeighbors(trainSetFish15,testSet1Fish15,personlist15,personlist1);
resultF15(2) = nearestNeighbors(trainSetFish15,testSet2Fish15,personlist15,personlist2);
resultF15(3) = nearestNeighbors(trainSetFish15,testSet3Fish15,personlist15,personlist3);
resultF15(4) = nearestNeighbors(trainSetFish15,testSet4Fish15,personlist15,personlist4);
resultF15(5) = nearestNeighbors(trainSetFish15,testSet5Fish15,personlist15,personlist5);

end


%% eiganfaces function
function[ureduce,utot] = eigenfaces(subset)
d = 9; %30 %9
[m,n] = size(subset);
x = {};

% sum all the 2500x1 faces in the subset
mu = sum(subset,2);

% find the mean face of them all
uj = (1/n).*mu;

% subtract mean from each image in subset
for i=1:n
    x{i} = (subset(:,i)-uj);
end

% unfold cells into x matrix
xmat = cell2mat(x);
% take SVD of it
[U,S,V] = svd(xmat);
%reduce to the first d eigenvectors
ureduce = U(:,1:d); %first d eigenvectors cutting out first 3

% find the 2500x1 variance
variance = diag(S).^2;

%find areas of max variance
% I = find(variance<d,1)
% sortedVar = sort(variance,'descend');

%% This is for plotting
% for k=1:9
%     subplot(3,3,k)
%     imagesc(reshape(U(:,k),50,50))
%     colormap gray
% end

ureduce = normc(ureduce);
utot = uj;
end

function[conf] = nearestNeighbors(trainSet,testSet,trainPerson,testPerson)

results = [];

%vector of people
per = zeros(1,10);

% find how many times each person appears in set
for j=1:10
    person = find(testPerson==j);
    numPer(j) = numel(person);
end


for i = 1:numel(testSet)
    for j = 1:numel(trainSet)
        %find the differences between the testset and training sets
        d(j) = norm(testSet{i}-trainSet{j});
    end
    % sort the differences
    [B,I] = sort(d);
    % take the 2 smallest differences
    d1 = B(1);
    d2 = B(2);
    % calculate the ratio between the two
    ratio = d1/d2;
    if ratio <1
        %if match, determine if people are the same
        if testPerson(i) == trainPerson(I(1))
            % if it is the same, add one to the correct person
            per(testPerson(i)) = per(testPerson(i))+1;
        end
    end
end

% calculate the percentages of the results based on how many times 
%each person showed up
for k = 1:10
    results(k) = 1-(per(k)/numPer(k));
end

% sum the results to get total correct percentage
conf = 1/10*sum(results);

end


%% fisherfaces algorithm
function[Wopt,utot] = fisherFaces(subset,personList)
%% PCA
c = 10; %10 and 31
k = c-1;
[mt,nt] = size(subset);
% sum all the 2500x1 faces in the subset
mu = sum(subset,2);

% find the mean face of subset
utot = (1/nt).*mu;

% subtract mean from each image in subset
% divide by standard deviation
for i=1:nt
    S = std(subset(:,i));
    x{i} = 1/S*(subset(:,i)-utot);
end

xmat = cell2mat(x);

%calculate the pca
[U,S,V] = svd(xmat);
% find the 2500x1 variance
variance = diag(S).^2;

%reduce to the first k eigenvectors
ureduce = U(:,1:k); %first k eigenvectors

% apply to unfolded x matrix
pcaMat = ureduce'*subset;
% now have eigenspace to apply FLD to.

%% average of eigenspace
[m,n] = size(pcaMat);
% sum all the reduced faces
mu = sum(pcaMat,2);
% find the mean face of them all
u = (1/m)*mu;

%% average of class data
C = 10; %10 classes
normface = [];
xprime = [];
Sdum = {};
for i=1:C
    %find all images of person i
    x = find(personList == i);
    n = numel(x); %number of times person i appears in subset
    for j=1:numel(x)
        % collect all faces of person i into an array
        % 2500 x how many images of person i
        normface = horzcat(normface,pcaMat(:,x(j)));
    end
    %    sum all the faces of person i
    mu = sum(normface,2);
    %    find the mean face of person i
    ui = (1/n)*mu;
    
    for y=1:n
        % subtract mean from normalized faces
        normalizedCells{y} = normc(normface(:,y)-ui);
    end
    
    % calculate scatter for class i
    for z=1:n
        Sdum{end+1} = (normalizedCells{z}*normalizedCells{z}');
    end
    
    %calculate Si
    Si{i}=sum(cat(3,Sdum{:}),3);
    
    %calculate xprime, not necessarily needed
    for element = 1:numel(normalizedCells)
        xprime{end+1} = normalizedCells{element};
    end
    
    %calculate Sb for class i
    SbCell{i} =  n*((ui-u)*(ui-u)');
end
% calculate Sb for class
Sb = sum(cat(3,SbCell{:}),3);

% calculate Sw for class
Sw = sum(cat(3,Si{:}),3);

% calculate Wfld through solving eigenvalue problem
[V,D,Wfld] = eig(Sb,Sw);

% calculate wopt'
Wopt = Wfld'*ureduce';

%untranspose wopt
Wopt = Wopt';

end

function[w] = applyUreduce(ureduce,subset,utot)
%% preprocessing
[m,n] = size(subset);
x = {};
for i=1:n
    S = std(subset(:,i));
    %subtract training mean from test subset
    x{i} = 1/S*(subset(:,i)-utot);
end
% unfold cells into x matrix
xmat = cell2mat(x);

%% apply ureduce to preprocessed subset
wmat = ureduce'*xmat;
w = num2cell(wmat,1);
end
