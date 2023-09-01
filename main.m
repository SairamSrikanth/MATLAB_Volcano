clear
num_scenario = 100;
[G_idx, rho, mu_idx, rc, M_idx] = randomVolcano(num_scenario);
path="VolcanoGen_Dataset/Volcano_10";
mkdir(path);

sigma=[0.1,0.25,0.5];
counter=0;
for i = 1:100
    synthetic_data = synthetic_data_generator_v3(G_idx, rho, mu_idx, rc, M_idx(i));
    if (synthetic_data(4, end)<0)
        disp("skipped");
        continue;
    end
    if (-synthetic_data(8,end)*10^9<1)
        disp("skipped");
        continue;
    end    
    for j = 1:3
        t = synthetic_data(1,:);
        tilt = synthetic_data(8,:);
        interval=randi([1000, 1500]);
        if length(t) > interval
            interval=length(t);
        end
        r = normrnd(0, sigma(j), [1, interval])*10^-9;
        padding = linspace(-interval+length(t)+1, 0,interval-length(t));
        t = [padding, t];
        t = t - t(end);
        padding = zeros(1, interval-length(tilt));
        tilt=[padding, -tilt];
        tilt=tilt+r;
        results=[t; tilt];
        plot(t, tilt);
        clf;
        counter=counter+1;
        writeData(path+"/observation"+counter, results, G_idx, rho, mu_idx, rc, M_idx(i), sigma(j), -synthetic_data(8, end)*10^9)
    end 
end