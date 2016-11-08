function tuss()

% Tuss the Truss simulator
% EK301
% Summer 2016
% Section C1
% Group : Anand, Konstantino, Vijay
% 27th July, 2016
%
% Analyses the truss design that the user inputs.
% Truss is expected to be in J= M - 3 configuraiton
% where j is the numebr of pin joints and M is the number of members.
% Three reaction forces are also considered at any joints. One has to be a
% pin joint and other a roller joint.
%
% Created by:
% Vijay H. Thakkar -- <(") (BU ENG 2019)
% 
% Summer 2016
% v.1.0
% This was a freaking pain in the ass to code. Least you could do it pass 
% it on to prevent others to go thought the same pain as it did :)
% Remember to chage the parameters a and b in section 5 to your specific lab
%
% ALL CAPS COMMENTS ARE NOTES FOR NEXT VERSION UPDATE


%% 1 - DEFINE VARIABLES (readme.txt)
clear
clc
% r is the numebr of joints and the rows in the A matrix
% c is the number of joints and the colums in the A matrix
% A is the matrix of r vs c and shows the connection between the joints and 
%     the truss members
% C is the porgram specific matrix that displays 1 for a connecton between
%     joint and member at that row cross column position correspoiding to A.
%     The purpose of this matrix is to ask for distances between the joints
%     that serve as the truss member lengths.
% MemLength is a row vector of c columns that stores the length of each
%     member indexed according to member number
% X is a vector that stores the X coordinates of the joints
% Y is a vector that stores the Y coordinates of the joints
% L is a colum vector that indexes the load at joint index in newtons
% a is the constant of multiplication for the semi-empirical fit
% b is the power to which straw length is raised for the semi-empirical fit
% MemLengths stores the calculated lenghts of each member
% cost stores the cost of the truss once calculated
% terminate is a logical true or false that stores the status of script
%     termination condition
% pinjoint is the number of the joiny which takes both x and y loads
% rollerjoint is the number of the number of the joint which takes only y
%     loads

%% 2 - Get truss size and and initiate variables

% NEED TO ADD A DISCLAIMER HERE AND INSTRUCTIONS ON HOW TO USE THE PORGRAM


%ask for truss size
fprintf('Please enter the number of joints in your truss below \n')
r = GetPosIntegerParameter();
c = 2*r-3;

%show input and ask if AOK
fprintf('The number of joints you have entered is %d \n',r)
fprintf('The number of truss memebrs is therefore M = 2J - 3 = %d \n',c)
fprintf('Are these parameters okay? Y/N --> ')

%get response
    str = GetYN();

%If the user types it worng they can change it
   while (str == 'N') || (str == 'n')
        clc
        r = GetPosIntegerParameter();
        c = 2*r-3;
        fprintf('The number of joints you have entered is %d \n',r)
        fprintf('The number of truss memebrs is therefore M = 2J - 3 = %d \n',c)
        fprintf('Are these parameters okay? Y/N --> ')
        str = GetYN();
   end
   clear str
clc      
%Otherwise initialize the  matricies and then  carry on

A = zeros(2*r);
C = zeros(r,2*r-3);
MemLength = zeros(1,c);
X = zeros(1,r);
Y = zeros(1,r);
L = zeros(2*r,1);
T = zeros(2*r,1);
cost = 0;
terminate = false;
pinjoint = 0;
rollerjoint = 0;

%% 3 - Get all the connections (matrix C)

fprintf('Your truss has %d joints and %d members. \n',r ,c)
fprintf('Now we will ask you which members are connected to which joints\n')

for i = 1:r
    fprintf('How many members connect to joint no. %d ?\n',i)
    % NoCons is a temp parameter within the for loop storing no.
    % of members connect to the ith joint. 
    NoCons = GetPosIntegerParameter();
    fprintf('Joint %d connects to %d truss members. Is this correct? Y/N -->',i, NoCons)
    str = GetYN();
    while (str == 'N') || (str == 'n')
          clc
          fprintf('\n')
          fprintf('Re-enter number of connected members here -->')
          NoCons = GetPosIntegerParameter();
          fprintf('Joint %d connects to %d truss members. Is this correct? --> Y/N',i, NoCons)
          str = GetYN();
    end
    clc
    clear str
    fprintf('Please enter below one by one\n')
    fprintf('Which members does joint no. %d connect to? \n',i)
    Cons = GetNInputs(i,NoCons);
    for j = 1:NoCons
        C(i,Cons(j))= true;
    end 
    clc
end
%Now draw the C matrix for the user
C

%Now we check that every colum has 2 ones in it. No more and no less
%Otherwise the parameters have been entered icorrectly

cons = sum(C);
if ~isempty(find(cons > 2))
    sentence = ['The connections you have entered do not agree with the premise' ...
    'of the truss desgin. Terminating program.'];
    error(sentence)
else
end

%% 4 - Get lenghts of the truss members
% pause(1)
% fprintf('Now we will ask for the lenghts of the members below\n')
% fprintf('Please enter the lengths in centemeters one by one\n')
% pause(1)
% %This for loop now gets the lenghts of all the members and stores them in
% %MemLenght indexed according to member number
% for i = 1:c
%     fprintf('Please enter the length of member %d. ',i)
%     TempLength = GetRealParameter();
%     fprintf('The length of member %d is %.4f cm\n',i ,TempLength)
%     fprintf('Is this the correct length? Y/N\n')
%     str = GetYN();
%     while (str == 'N') || (str == 'n')
%           clc
%           TempLength = GetRealParameter();
%           fprintf('The length of member %d is %.4f cm. Ok? Y/N -->\n',c ,TempLength)
%           str = GetYN();
%     end
%     clear str
%     %Now devide length by 100 to convert to meters for the HOLY SI UNIT
%     %STANDARD!!! FUCK  IMPERIAL UNITS!!!!
%     MemLength(i) = TempLength/100;
%     clc
%     MemLength
% end

%% 5 - Get the reference Coordinates of joints

%X coordinates first
pause(1)
fprintf('Now we will ask for the reference coordinates of all the joints\n')
fprintf('Please enter the coordinates in centemeters one by one\n')
fprintf('First we will get the X coordinates \n')
for i = 1:r
    %TempCor stores the temporary coordinate before confirmation
    fprintf('\n')
    fprintf('Please enter the X coordinate of joint %d. ',i)
    TempCor = GetRealParameter();
    fprintf('X Coordinate of joint %d is %.4f. Is this correct? Y/N? \n',i,TempCor)
    str = GetYN();
    
    while (str == 'N') || (str == 'n')
        clc
        fprintf('Please enter the X coordinate of joint %d. ',i)
        TempCor = GetRealParameter();
        fprintf('X Coordinate of joint %d is %.4f. Is this correct? Y/N? \n',i,TempCor)
        str = GetYN();
    end
    
    X(i) = TempCor;
    clear str
    clc
end
pause(1)

%Now for Y coordinates
fprintf('Now we will ask for the Y corodinates. Please enter one by one \n')
pause(1)

for i = 1:r
    %TempCor stores the temporary coordinate before confirmation
    fprintf('Please enter the Y coordinate of joint %d. ',i)
    TempCor = GetRealParameter();
    fprintf('Y Coordinate of joint %d is %.4f. Is this correct? Y/N? \n',i,TempCor)
    str = GetYN();
    while (str == 'N') || (str == 'n')
        clc
        fprintf(' \n Please re-enter the Y coordinae of joint %d --> ',i)
        TempCor = GetRealParameter();
        fprintf(' \n X Coordinate of joint %d is %.4f. Is this correct? Y/N? -->',i,TempCor)
        str = GetYN();
    end
    Y(i) = TempCor;
    clear str
    clc
end

%% 6 - Get loaded joint and value of load

fprintf('Now we will ask for the loaded joints number. \n')
fprintf('Please enter the number of the joint at which the truss is loaded. \n')
pause(1)

%First get the joint at which the truss is loaded
LoadJoint = GetPosIntegerParameter();
fprintf('The load is at joint %d Is this Correct? Y/N --> \n',LoadJoint)
str = GetYN();
while (str == 'N') || (str == 'n')
        clc
        fprintf(' \n Please re-enter the number of the joint at which the truss is loaded --> ')
        LoadJoint = GetPosIntegerParameter();
        fprintf(' \n The load is at joint %d Is this Correct? Y/N --> ',LoadJoint)
        str = GetYN();
end
clear str

%Now get the value of load in Newtons
fprintf('Now enter the load at joint %d. \n',LoadJoint)
Load = GetRealParameter();
fprintf('The load at joint %d is %.4f N. Is this Correct? Y/N --> \n',LoadJoint, Load)
str = GetYN();
while (str == 'N') || (str == 'n')
        clc
        fprintf('Please re-enter the load at joint %d ',LoadJoint)
        Load = GetRealParameter();
        fprintf('\n The load at joint %d is %.4f N. Is this Correct? Y/N --> ',LoadJoint,Load)
        str = GetYN();
end
clear str
clc
%Index load in the L matrix in 2*joint column to represnt vertical load and
%then mmultiply by -1 to represent load in -ve direction
xcor = LoadJoint + r;
L(xcor,1) = Load*-1;
clear xcor
pause(1)

%% 7 - Get 2 reaction force joint numbers
fprintf('Now enter the joints which bare the reaction forces \n')
%Pin Joint First
fprintf('First enter location of the pin joint here. ')
    pinjoint = GetPosIntegerParameter();
    fprintf('Pin Joint is joint no. %d. Is this correct? Y/N -->  ',pinjoint)
    str = GetYN();
    while (str == 'N') || (str == 'n')
        clc
        fprintf('Please re-enter the pin joint no. here. ')
        pinjoint = GetPosIntegerParameter();
        fprintf('Pin Joint is joint no. %d. Is this correct now? Y/N --> ',pinjoint)
        str = GetYN();
    end
    clc
    clear str

%Roller Joint now
pause(1)
fprintf('Now enter location of the roller joint here. ')
    rollerjoint = GetPosIntegerParameter();
    fprintf('Roller joint is joint no. %d. Is this correct? --> ',rollerjoint)
    str = GetYN();
    while (str == 'N') || (str == 'n')
        clc
        fprintf('Please re-enter the roller joint no. here --> ')
        rollerjoint = GetPosIntegerParameter();
        fprintf('Roller Joint is joint no. %d. Is this correct now? Y/N --> ',rollerjoint)
        str = GetYN();
    end
    clear str
    
%Check is the pin joint isnt the same as the roller joint
if pinjoint == rollerjoint
    error('Pin joint cannot be the same as the roller joint')
end

pause(1)

%% 8 - Get Regression Parameters

% NEED TO UNBAKE THE A AND B PARAMETER HERE AS WELL AS IN THE LATER SECTIOINS

% fprintf('Now we will ask for the regression parameters of the load tests \n')
% fprintf('Please enter one by one below. Note that all loads are in Newtons \n')
% fprintf('And that all lengths of memebers are in centimeters \n')
% fprintf('Type of fit is a semi empirical fit of the data form lab 1 \n')
% fprintf('Fitting equation is of the form F(x) = a * x.^(-b) \n')

% fitting parameters for EK301 C1, summer 2016= 
% 
%      General model:
%      f(x) = (a*x.^b)
%      Coefficients (with 95% confidence bounds):
%        a =         236  (146.8, 325.2)
%        b =       -1.266  (1.104, 1.428)

% %First get slope of line
% fprintf('First enter the value of a here. ')
%     a = GetRealParameter();
%     fprintf('\n Value of a is %.4f. Is this correct? Y/N --> ', a)
%     str = GetYN();
%     while (str == 'N') || (str == 'n')
%         clc
%         fprintf('\n Please re-enter the slope value here. ')
%         a = GetRealParameter();
%         fprintf('\n Slope value is %.4f. Is this correct? Y/N --> ', a)
%         str = GetYN();
%     end
%     clear str
% pause(1)

% %Now get b
% fprintf('Now enter the value of b here. ')
%     b = GetRealParameter();
%     fprintf('\n Value of b is %.4f. Is this correct? Y/N --> ', b)
%     str = GetYN();
%     while (str == 'N') || (str == 'n')
%         clc
%         fprintf('\n Please re-enter the intercept value here --> ')
%         b = GetRealParameter();
%         fprintf('\n intercept value is %.4f. Is this correct? Y/N --> ', b)
%         str = GetYN();
%     end
% clc
% clear str
% pause(1)

a = 236;
b = -1.266;

%% 9 - Now do all the calculations!!! -- first A and then T

%First calculate lengths of all members based on joint numbers
MemLengths = CalcLengths(X, Y, C);

fprintf('Now we will calcualte the cost of the entire truss \n')
[cost, terminate] = CalcCost(r, MemLengths);
if terminate == true
    error('lel u suk m8. k bye')
elseif terminate == false
    fprintf('Okay we will proceed with the load analysis now. \n')
    fprintf ('Please wait for the results to be calculated. \n')
end

% First make the A matrix
% Ill have you know, writing the code for making A took about an hour itself
% following an intsnt debugging session of another hour
% I suggest trying to code it yourself

% Who am i kidding here. You are now. There is a reason why you have my code
% with you lol.

A = MakeMatrixA(C, X, Y, MemLengths, pinjoint, rollerjoint);

T = CalculateForces(A, L);

%Pass on T to a funtion that displays all forces in members in a near
%setence structure and indicated weather member is in tension or compression

DisplayTensions(T)

[delta, failstate , Maxload] = CalculateDelta(T, MemLengths, a, b);

DisplayDelta(delta)

T = T';
Tloads = T(1:(length(T)-3));
Treactions = T((length(T)-2):length(T));
CompressionIndex = find(Tloads > 0);
Coms = Tloads(CompressionIndex);
FailedCom = max(Coms);
FailedMemtemp = find(Tloads == FailedCom);

if length(FailedMemtemp) > 1
    FailedMem = FailedMemtemp(1);
else
    FailedMem = FailedMemtemp;
end

MemLengthOfFailedBeam = MemLengths(FailedMem);
MaxLoadItCouldTake = Maxload(FailedMem);
UpMax = 325.2*MemLengthOfFailedBeam^(-1.266);
PerUpMax=(100*(UpMax-MaxLoadItCouldTake)/MaxLoadItCouldTake);
ratio = (max(Coms)/Maxload(FailedMem));

if failstate == 1 
    fprintf('Truss does not survive the load of %.4f N. \n', Load)
    fprintf('Member %d buckles first by a delta load of %.4f N \n', FailedMem, delta(FailedMem))
    fprintf('Ratio of actual load to max compressive load is %.4f \n', ratio)
    fprintf('Maximum compression that hte member can take is %.4f N \n', MaxLoadItCouldTake);
    fprintf('Uncertainty in this buckling load prediction is %.2f percent \n', PerUpMax);
elseif failstate == 0
    fprintf('Truss survives the load of %.4f N.\n', Load)
    fprintf('Ratio of actual load to max compressive load is %.4f \n', ratio)
    fprintf('Member %d would buckle first \n', FailedMem)
    fprintf('Maximum compression that the member can take is %.4f N \n', MaxLoadItCouldTake);
    fprintf('Uncertainty in this buckling load prediction is %.2f percent \n', PerUpMax);
end  


MaxTrussLoad = (Load/ratio);
UncMaxTrussLoad = MaxTrussLoad*PerUpMax/100;
CostToLoad = cost/MaxTrussLoad;
fprintf('Max load on truss is therefore %.4f N \n', MaxTrussLoad )
fprintf('Uncertainty in this max load on truss is %.2f N \n', UncMaxTrussLoad)
fprintf('Cost to load ratio is therefore %.2f Dolars per Newton \n', CostToLoad)
    %Now ask user if they want it in an excel file


Tmax = T/ratio;

fprintf('Save data in an excel file? Y/N --> \n')
str = GetYN();
if str == 'y' || str == 'Y'
    name = input('Enter name of file here: ','s');
    name = [name '.xlsx'];
    concat = {Tloads, Treactions, X, Y, MemLengths, Maxload, delta, cost, MaxTrussLoad, UncMaxTrussLoad,Tmax, CostToLoad};
    Ncols = max(cellfun(@numel,concat));
    array = nan(numel(concat),Ncols);
    for iRol = 1:numel(concat),
        array(iRol,1:numel(concat{iRol})) = concat{iRol};
    end
    xlswrite(name, array, 1)
end


end

%% 10 - Auxilliary Functions

% NEED TO MAKE THE FUNCTIONS IMMUNE TO EMPT INPUTS
% NEED TO MAKE A FUNCTION FOR ARRAY INPUTS 
% NEED TO CORRECT THE POS/NED SIGNS OF LOADS CALCULATED\


function out = GetPosIntegerParameter()
%Using b as a local substiture for r to emphesize scope differrence
b = input('Enter the number here --> ');
    while abs(floor(b))~=b
         b = input('Positive intigers only please. Re-eneter here --> ');
    end
out = b;
end

function out = GetRealParameter()
    out = double(0);
    out = input('Enter your number here --> ');  
end

function str = GetYN()
    str = input('','s');
    if isempty(str)
        str = 'y';
    else
        while ((str ~= 'y') && (str ~= 'Y') && (str ~= 'n') && (str ~= 'N'))
           str = input('Y/N responses only. Re-enter here --> ','s');
        end
    end
end

function varargout = GetNInputs(JNum,NoCons)
%Get a varibale (NoCons) number of inputs from the user and outs them in the
% fist cell of the struct called varargout
varargout{1} = zeros(1,NoCons);
    for i = 1:NoCons
        %Get the connection
        TempCon = GetPosIntegerParameter();
        %Confirm Connection
        fprintf('You said joint %d is connected to member %d.\n', JNum, TempCon)
        fprintf('Are you sure? Y/N -->')
        str = GetYN();
            while (str == 'N') || (str == 'n')
                TempCon = GetPosIntegerParameter();
                fprintf('You said joint %d is connected to member %d.',...
                'Are you sure \n',i,TempCon)
                str = GetYN();
            end
        varargout{1}(i) = TempCon;
        clear str
    end
end

function [cost, terminate] = CalcCost(joints, MemLengths)

TotalLength = sum(MemLengths);
% Calculates the cost of the truss and displays if cost is feasable and how
% much under or over budget the truss is
cost = (10*joints) + (1*TotalLength);
delta = abs(335-cost);

%overbudget process
if cost > 335
    fprintf('The cost of the truss is $ %.4f which is $ %.4f over budget \n', cost, delta)
    fprintf('Proceed with the load analysis? Y/N --> ')
    str = GetYN();
    if (str == 'N') || (str == 'n')
        terminate = true;
    else
        terminate = false;
    end
end

%underbudget process
if cost <= 335
    fprintf(' The cost of the truss is $ %.4f which is $ %.4f under budget \n', cost, delta)
    fprintf('Proceed with the load analysis? Y/N --> ')
    str = GetYN();
    if (str == 'N') || (str == 'n')
        terminate = true;
    elseif (str == 'Y') || (str == 'y')
        terminate = false;
    end
end

end

function A = MakeMatrixA(C, X, Y, MemLengths, pinjoint, rollerjoint)
%Makes the A matrix that needs to be inverted for load calculation
%Setup parameters for the loopoopoopooopoopiopopopopopopoppoo
[joints, ~] = size(C);
members = joints*2-3;
r = joints;
A = zeros(r*2);
%index reaction force vectors first

A(pinjoint,2*r-2) = 1;
A((pinjoint+r),2*r-1) = 1;
A((rollerjoint+r), 2*r) = 1;

%index member force vectors now if the joints are connected
for i=1:r
    for j=1:members
        if C(i,j) == true
           tempvec = find(C(:,j) == true);
           tempindex = find(tempvec ~= i);
           oj = tempvec(tempindex);
           xcor = i;
           ycor = i + r;
           A(xcor,j) = (X(oj)-X(i))/MemLengths(j);
           %Rinse and repeat for y direction
           A(ycor,j) = (Y(oj)-Y(i))/MemLengths(j);
        elseif C(i,j) == false
           %if not connected then make the entry 0
           xcor = i;
           ycor = i + r;
           A(xcor,j) = 0;
           A(ycor,j) = 0;
        end
    end
end
end

function T = CalculateForces(A, L)
%This calcualtes the forces on each member
T = (inv(A))*(L);
end

function DisplayTensions(T)
%This funtion displays all calculated forces in a neat manner and shows if
%memebr is under tension or compression

for i = 1:(length(T)-3)
    if T(i) > 0
       %Compressive load case
       fprintf('Member number %d is under a compressive load of %.4f Newtons. \n',i ,abs(T(i)))
    elseif T(i) < 0
       %Tensile load case
       fprintf('Member number %d is under a tensile load of %.4f Newtons. \n',i ,abs(T(i)))
    elseif T(i) == 0
       %No load at all case
       fprintf('Member number %d is under no load at all. \n',i )
    end
end

for i = 1:3
    fprintf('Reaction force S %d is %.4f Newtons\n',i ,abs(T(length(T)-3+i)))
end
    fprintf('S1 is the X directional load at pin joint. \n')
    fprintf('S2 and S3 are the Y directional load at pin joint and roller joint respectively. \n')
end

function [Delta, failstate, MaxLoad] = CalculateDelta(T, MemLengths, a, b)
%This function calcuatles what the delta load is between the predicted max
%compressive for a member can take based on the regression fit curve and
%the actual load under the specified weight on truss
%It also retuns if the truss breaks or not and if yes at which memebr first

MaxLoad = a*(MemLengths).^b;

%we only have to  test for compression
%Failstate is a logical variable that stores true if the truss fails for the
%given load nad flase if it doesnt.
T = T';
Tloads = T(1:(length(T)-3));
%if truss fails lol 
CompressionIndex = find(Tloads > 0);
Compressions = abs( Tloads(CompressionIndex) );
MaxCompressions = abs( MaxLoad(CompressionIndex) );

if sum(Compressions > MaxCompressions) == 0
    failstate = false;
elseif sum(Compressions > MaxCompressions) > 0
    failstate = true;
end

%If it doesnt fail lol
if failstate == false
    fprintf('Congratulations! The truss survives! \n')
    Delta = abs(Tloads) - abs(MaxLoad);
elseif failstate == true
    fprintf('The truss fails with the given load. Try a different value. \n')
    Delta = abs(MaxLoad) - abs(Tloads);
end

end

function DisplayDelta(delta)
%This funtion displays all delta values between max loads and actual load
%on each memebr

for i = 1:length(delta)
    fprintf('Delta load for compression member no. %d is %.4f newtons. \n',i ,delta(i))
end

end

function MemLengths = CalcLengths(X, Y, C)

mems = length(X)*2-3;
MemLengths = zeros(1,mems);

for i = 1:mems
    tempindex = find(C(:,i) == true);
    j1 = tempindex(1);
    j2 = tempindex(2);
    MemLengths(i) = sqrt(((X(j1)-X(j2))^2)+((Y(j1)-Y(j2))^2));
    clear tempindex tempvec
end
end
