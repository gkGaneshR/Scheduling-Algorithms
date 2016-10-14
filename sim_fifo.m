function out=sim_fifo()
% Simulation
% The FIFO Scheduling Algorithm

% switch size
N = 16;

% traffic load
rho = [0.2:0.02:0.7];

% average delay
Avg_delay_fifo = zeros(1,length(rho));

% FIFO queue size
fifo_len = 1000;

% total simulation time
T = 50000;

for r = 1:length(rho);
    
    % start discrete time random event simulation
    t = 1;
    
    % FIFO queue length
    fifo_l = zeros(1,N);    
    
    % input FIFO queue
    input_fifo = zeros(N,fifo_len);
    
    % FIFO queue packet arrival time
    fifo_arr = zeros(N,fifo_len);
    
    % total latency
    total_delay = 0;
    % total transmitted packets
    total_pkt = 0;
    
    while t<=T
        
        % packets arrive at inputs (if pk_arr(i)=1, then there is a packet 
        % at input i; if pk_arr(i)=0, then no packet at input i)
        pk_arr = rand(1,N)<=rho(r);
        
        % destined output ports
        output_pt = ceil(N*rand(1,N)).*pk_arr;
        
        % push packets into input queues
        for j=1:N   % for every input j
            
            if output_pt(j)>0  % if input j has a packet arrival
                % ------FIFO
                input_fifo(j,fifo_l(j)+1) = output_pt(j);
                
                % store the arrival time
                fifo_arr(j,fifo_l(j)+1) = t;
                
                % increase the queue length
                fifo_l(j) = fifo_l(j) + 1;              
                                            
            end
            
        end        
        
        % FIFO scheduling -----
        % HOL packets
        HOL_fifo = input_fifo(:,1);
        
        for u=1:N % check if there are packets destined to output u
            if sum(HOL_fifo==u)==1  % only one packet going to one output
                
                % temp_i returns the input that has a packet to output u
                temp_i = find(HOL_fifo==u);
                
                % shift the queue forward by one packet
                % delete the transmitted packet
                temp_q = input_fifo(temp_i,[2:end]);                
                input_fifo(temp_i,[1:end-1]) = temp_q;                
                               
                % calculate delay
                total_delay = total_delay + (t-fifo_arr(temp_i,1));
                
                % delete the arrival time of the transmitted packet
                temp_p = fifo_arr(temp_i,[2:end]);                
                fifo_arr(temp_i,[1:end-1]) = temp_p;
                
                % FIFO queue length decreases by 1
                fifo_l(temp_i) = fifo_l(temp_i) - 1;
                
                % calculate total packets transmitted
                total_pkt = total_pkt + 1;
                
            elseif sum(HOL_fifo==u)>1  % more than one packet going to one output
                
                % temp_k returns those inputs who have packets to output u
                temp_k = find(HOL_fifo==u);
                
                % randomly select an input to send output u
                rnd_k = rand(1,length(temp_k)); % generate random numbers                
                temp_i = temp_k(find(rnd_k==max(rnd_k))); % pick the one with the largest random number
                
                % shift the queue forward by one packet
                % delete the transmitted packet
                temp_q = input_fifo(temp_i,[2:end]);                
                input_fifo(temp_i,[1:end-1]) = temp_q;
                
                % calculate delay
                total_delay = total_delay + (t-fifo_arr(temp_i,1));
                
                % delete the arrival time of the transmitted packet
                temp_p = fifo_arr(temp_i,[2:end]);                
                fifo_arr(temp_i,[1:end-1]) = temp_p;
                
                % FIFO queue length decreases by 1
                fifo_l(temp_i) = fifo_l(temp_i) - 1;
                
                % calculate total packets transmitted
                total_pkt = total_pkt + 1;                
                
            end
        end  
        
        t = t+1;        
    end
    
    % average cell latency
    Avg_delay_fifo(r) = total_delay/total_pkt;
    
end

out=[rho;Avg_delay_fifo];
% save fifo.mat out;
% 
% figure(1)
% semilogy(rho,Avg_delay_fifo);
% axis([0.2 1 1e-1 1e+3]);

