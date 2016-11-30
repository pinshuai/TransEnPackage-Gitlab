%RunMetabolism.m
%This is the run code for a Shannon entropy analysis of a metabolism
%dataset. It is configured to specific columns of the input file, and it
%performs some additional data pre-processing. However, it should be readily
%adaptable by the user to analyze other inputs of interest. The run code
%should determine which pairwise combinations of variables should be run
%in the mutual information and transfer entropy functions, as well as which
%outputs should be saved.

%USER INPUT SECTION
A = importdata('/Users/lglarsen/Dropbox/LaurelsTransEnToolbox/Data/AccDailyAll.csv', ',');%Enter the complete path to the spreadsheet containing the metabolism data.
A.textdata{26} = 'Discharge, cfs, range';
alpha = 0.05;% Significance level
nbins = [11, 11, 11]; %Nubmer of bins to use in 1, 2, and 3 dimensions

%PERFORM ADDITIONAL DATA PRE-PROCESSING
A.data(:,26) = A.data(:,23)-A.data(:,24); %Range in discharge
get_rid_of_bad_conductivity = A.data(:,10:13);
get_rid_of_bad_conductivity(get_rid_of_bad_conductivity<5) = NaN; %Get rid of datapoints lower than DI!
A.data(:,10:13) = get_rid_of_bad_conductivity;

%INITIALIZE MATRICES
Imat = NaN(3, 24); %One extra one: discharge range
Icritmat = Imat;
Tfirstmat = Imat;
Tcritmat = Imat;
Tbiggestmat = Imat;
Ishortmat = Imat;
Icritshortmat = Imat;

%LOOP OVER ALL PAIRS OF SOURCE AND SINK VARIABLES TO CALCULATE MUTUAL
%INFORMATION AND TRANSFER ENTROPY
for ii = 3:5 %Sink nodes (information receivers)
    for jj = 3:26 %Source nodes (information transferrers)
        if jj ~= ii %Don't do it if source = sink
            M = A.data(:,[jj, ii]); %Selectthe pair of variables to run
            if jj<=5, M(:,1) = sqrt(abs(M(:,1))); end %Perform transformation to an approximately normal distribution
            if ii<=5, M(:,2) = sqrt(abs(M(:,2))); end %Perform transformation to an approximately normal distribution
            if jj>=10, M(:,1) = log(abs(M(:,1))); end %Perform transformation to an approximately normal distribution for non-temperature source variables.
            sink = A.textdata{ii}; %Text descriptor of sink node
            source = A.textdata{jj};%Text descriptor of source node
            I = mutinfo(M, nbins); %Calculate mutual information to save
            Imat(ii-2, jj-2) = I; %Save it in a matrix
            Icrit= mutinfo_crit(M, nbins, alpha); %Calculate the significance threshold for mutual information
            Icritmat(ii-2, jj-2) = Icrit; %Save it in a matrix
            T = NaN(1, 181); %Initialize the vector of transfer entropy over the range of lags examined (here, 180 days)
            Tcrit = T; %Initialize the vector of the transfer entropy significance threshold over the range of lags examined.
            for lag = 0:180 %Run all lags that will still generate 500 overlapping datapoints (if no blanks). Here, we are running lags up to 180 days.
                [t, N, Mshort] = transen(M, lag, nbins); %Calculate transfer entropy from the input data vectors over the range of lags up to and including "lag"
                if N >=500 %Minimum number of valid samples for saving transfer entropy
                    T(lag+1) = t; %Save the calculated transfer entropy
                    Tcrit(lag+1) = transen_crit(Mshort, lag, alpha, 1000, nbins); % Compute the significance threshold and save it
                end
            end
            if length(find(T>=Tcrit)) >=1, Tfirstmat(ii-2,jj-2) = T(find(T>=Tcrit, 1)); end %Save the biggest value of T over the significance threshold
            Tbiggestmat(ii-2,jj-2) = T(find(T-Tcrit==max(T-Tcrit),1)); %Save the first value of T over the significance threshold
            figure
            plot(0:180, T, 'b-'), xlabel('Lag, days'), ylabel('Tz')
            hold on 
            plot(0:180, Tcrit, 'b--')
            title(sprintf('%s%s%s', source, '-->', sink)), pause(0.1)
            savefig(sprintf('%s%i%s%i', '/Users/lglarsen/Dropbox/LaurelsTransEnToolbox/Results/Accotink_t-1untransformed11bins', jj, ' to ', ii)) %Save the output graphic
        end
    end
end
            
