clear all
close all
clc
adcData=readDCA1000('adc_1GHZ_position3_ (1).bin',200,4);
numADCSamples=342;
numChirps=900;
virsual_Rx=12;

data=zeros(numChirps,numADCSamples,virsual_Rx);
for i=1:virsual_Rx
    data(:,:,i)=reshape(adcData(i,:) , numADCSamples,numChirps ).';
    % data(:,:,i)=data(:,:,i).';
end
data(:,:,5:8)=[];%Obtain the Tx1 and Tx3's data only.
save two_people_1.mat data; %Save the data named two_people_1.mat.
