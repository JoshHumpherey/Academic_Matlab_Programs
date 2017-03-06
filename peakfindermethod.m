%% Trying to find freq peaks with findpeak function

%% initializing the variables
i = 0;
j = 0;

%%inverting the data to utilize the findpeak function
for j = 1:1250
    for i = 1:666
        inverted(i,j) =(-1)*columns(i,j);
    end
end

%%Smoothing the data to remove false positives
for j = 1:1250
    for i = 1:666
        smoothinverted(i,j) = smooth(inverted(i,j), 'moving');
    end
end

%%adding in predeterimined lengths
extra0 = zeros(4,1250);
extra1 = zeros(3,1250);
extra2 = zeros(2,1250);
extra3 = zeros(1,1250);


%% finding the peaks of the inverted data
for j = 1:1250
    final = (smoothinverted(:,j));
    z = length(findpeaks(final));
    extra = 4-z;

        if extra == 0
            extra0(:,j) = findpeaks(final);
        elseif extra == 1
            extra1(:,j) = findpeaks(final);
        elseif extra == 2
           extra2(:,j) = findpeaks(final);
        elseif extra == 3
            extra3(:,j) = findpeaks(final);
        end      
end

%%Resizing the matrices
for j = 1:1250
    z = length(extra1(:,j));
    if z < 4
        extra1(4,j) = 0;
    end
end

for j = 1:1250
    z = length(extra2(:,j));
    if z < 4
        extra2(3,j) = 0;
        extra2(4,j) = 0;
    end
end

for j = 1:1250
    z = length(extra3(:,j));
    if z < 4
        extra3(2,j) = 0;
        extra3(3,j) = 0;
        extra3(4,j) = 0;
    end
end

%%combining the resized matrices
FinalCombMin = extra0 + extra1+extra2+extra3;

%%replacing zeros with the parent frequency
for j = 1:1250
    if FinalCombMin(2,j) == 0
        FinalCombMin(2,j) = FinalCombMin(1,j);
    end
end


%%Re-inverting the data to compare with original set
FinalCombMin = -FinalCombMin;

%%Splitting Final Combination into two parts for easier freq manipulation
FinalFlatMin = zeros(4,1250);
for j = 1:1250
    FinalFlatMin(1,j) = FinalCombMin(1,j);
    FinalFlatMin(2,j) = FinalCombMin(2,j);
    FinalFlatMin(3,j) = FinalCombMin(3,j);
    FinalFlatMin(4,j) = FinalCombMin(4,j);
end

%%Relating the inverted S(1,1) values back to their peak frequency
row1min = zeros(1,1250);
row2min = zeros(1,1250);
row3min = zeros(1,1250);
row4min = zeros(1,1250);

for j = 1:1250
    for i = 1:666
        if FinalFlatMin(1,j) == moddata(i,j)
            row1min(j) = i;   
        end
    end
end

for j = 1:1250
    for i = 1:666
        if FinalFlatMin(2,j) == moddata(i,j)
            row2min(j) = i;
        end
    end
end

for j = 1:1250
    for i = 1:666
        if FinalFlatMin(3,j) == moddata(i,j)
            row3min(j) = i;   
        end
    end
end

for j = 1:1250
    for i = 1:666
        if FinalFlatMin(4,j) == moddata(i,j)
            row4min(j) = i;
        end
    end
end

%% Finding the peak values as vectors P1 and P2
peak1min = zeros(1,1250);
peak2min = zeros(1,1250);
peak3min = zeros(1,1250);
peak4min = zeros(1,1250);

peak1min = 1 + row1min*(0.003)-0.003;
peak2min = 1 + row2min*(0.003)-0.003;
peak3min = 1 + row3min*(0.003)-0.003;
peak4min = 1 + row4min*(0.003)-0.003;

%% Removing the false positives
for j = 1:1231
        if peak3min(j) == 0.9970
            peak3min(j) = 0;
        end
        if peak4min(j) == 0.9970
            peak4min(j) = 0;
        end
end

%% Finding the maximums for the data
SmoothNormal = -(smoothinverted);

%%adding in predeterimined lengths
extra0n = zeros(4,1250);
extra1n = zeros(3,1250);
extra2n = zeros(2,1250);
extra3n = zeros(1,1250);


%% finding the peaks of the inverted data
for j = 1:1250
    final = (SmoothNormal(:,j));
    z = length(findpeaks(final));
    extra = 4-z;

        if extra == 0
            extra0n(:,j) = findpeaks(final);
        elseif extra == 1
            extra1n(:,j) = findpeaks(final);
        elseif extra == 2
           extra2n(:,j) = findpeaks(final);
        elseif extra == 3
            extra3n(:,j) = findpeaks(final);
        end      
end

%%Resizing the matrices
for j = 1:1250
    z = length(extra1n(:,j));
    if z < 4
        extra1n(4,j) = 0;
    end
end

for j = 1:1250
    z = length(extra2n(:,j));
    if z < 4
        extra2n(3,j) = 0;
        extra2n(4,j) = 0;
    end
end

for j = 1:1250
    z = length(extra3n(:,j));
    if z < 4
        extra3n(2,j) = 0;
        extra3n(3,j) = 0;
        extra3n(4,j) = 0;
    end
end

%%combining the resized matrices
FinalCombMax = extra0n + extra1n +extra2n + extra3n;

%%replacing zeros with the parent frequency
for j = 1:1250
    if FinalCombMax(2,j) == 0
        FinalCombMax(2,j) = FinalCombMax(1,j);
    end
end

%%Splitting Final Combination into four parts for easier freq manipulation
FinalFlatMin = zeros(4,1250);
for j = 1:1250
    FinalFlatMax(1,j) = FinalCombMax(1,j);
    FinalFlatMax(2,j) = FinalCombMax(2,j);
    FinalFlatMax(3,j) = FinalCombMax(3,j);
    FinalFlatMax(4,j) = FinalCombMax(4,j);
end

%%Relating the inverted S(1,1) values back to their peak frequency
row1max = zeros(1,1250);
row2max = zeros(1,1250);
row3max = zeros(1,1250);
row4max = zeros(1,1250);

for j = 1:1250
    for i = 1:666
        if FinalFlatMax(1,j) == moddata(i,j)
            row1max(j) = i;   
        end
    end
end

for j = 1:1250
    for i = 1:666
        if FinalFlatMax(2,j) == moddata(i,j)
            row2max(j) = i;
        end
    end
end

for j = 1:1250
    for i = 1:666
        if FinalFlatMax(3,j) == moddata(i,j)
            row3max(j) = i;   
        end
    end
end

for j = 1:1250
    for i = 1:666
        if FinalFlatMax(4,j) == moddata(i,j)
            row4max(j) = i;
        end
    end
end

%% Finding the peak values as vectors P1 and P2
peak1max = zeros(1,1250);
peak2max = zeros(1,1250);
peak3max = zeros(1,1250);
peak4max = zeros(1,1250);

peak1max = 1 + row1max*(0.003)-0.003;
peak2max = 1 + row2max*(0.003)-0.003;
peak3max = 1 + row3max*(0.003)-0.003;
peak4max = 1 + row4max*(0.003)-0.003;

%%removing the false positives
for j = 1:1231
        if peak3max(j) == 0.9970
            peak3max(j) = 0;
        end
        if peak4max(j) == 0.9970
            peak4max(j) = 0;
        end
end

%% Placing all rows into two matrices
MinPeaks = zeros(4,1250);
MaxPeaks = zeros(4,1250);


MinPeaks(1,:) = peak1min;
MinPeaks(2,:) = peak2min;
MinPeaks(3,:) = peak3min;
MinPeaks(4,:) = peak4min;

MaxPeaks(1,:) = peak1max;
MaxPeaks(2,:) = peak2max;
MaxPeaks(3,:) = peak3max;
MaxPeaks(4,:) = peak4max;


figure
subplot(2,2,[1,2]);
plot(MinPeaks);
title('Minimum Peaks')

subplot(2,2,[3,4]);
plot(MaxPeaks);
title('Maximum Peaks')    

plot(MinPeaks);
title('Overlapping Peaks');
hold on
plot(MaxPeaks);
hold off
