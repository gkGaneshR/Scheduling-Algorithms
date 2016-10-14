function out=sim_rrm()
% Simulation
% The RRM Scheduling Algorithm

% switch size
N = 16;

% traffic load
rho = [0.2:0.02:0.7];

% average delay
Avg_delay_rrm = zeros(1,length(rho));

% VOQ queue size
voq_len = 100;

% total simulation time
T = 50000;

for r = 1:length(rho);
    
    % start discrete time random event simulation
    t = 1;
    
    % VOQ queue length
    voq_l = zeros(N,N);    
    
    % input VOQ queue
    input_voq = zeros(N,N,voq_len);
    
    % VOQ queue packet arrival time
    voq_arr = zeros(N,N,voq_len);
    
    % grant pointer for output
     gra_pointer = ones(1,N);
    
    % accept point for input
     acc_pointer = ones(1,N);    
    
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
        
        % push packets into input VOQ queues
        for j=1:N   % for every input j
            
            if output_pt(j)>0  % if input j has a packet arrival
                % VOQ
                input_voq(j,output_pt(j),voq_l(j,output_pt(j))+1) = output_pt(j);
                
                % store the arrival time
                voq_arr(j,output_pt(j),voq_l(j,output_pt(j))+1) = t;
                
                % increase the queue length
                voq_l(j,output_pt(j)) = voq_l(j,output_pt(j)) + 1;              
                                            
            end
            
        end        
        
        % RRM scheduling
        
        % Initialize request
        Req_rrm = zeros(N,N);
        
        % Initialize grant
        Gra_rrm = zeros(1,N);
        
        % generate requests
        for y = 1:N % for every input y, send requests to output
            
            % RRM Request: '1' means request, '0' means no request
            
             Req_rrm(y,:) = voq_l(y,:)>0;
           
        end
        
        % generate grants
        for u = 1:N  % for every output u, randomly select an input for grant
            if sum(Req_rrm(:,u))>0  % has at least one Request
                
                % select the input that is the closest to the grant pointer
                rnd_k = rand(1,N).*Req_rrm(:,u)'; % generate random numbers                
                temp_u = find(rnd_k==max(rnd_k)); % find the maximum random number
				
				% store the selected input in the Grant
                gra_pointer(u)=temp_u;
				
                % store the selected input into Gra_rrm(u)
                Gra_rrm(u) = gra_pointer(u);                
            end
        end
        
        % accept grants
        for y = 1:N   % for every input y, select a grant
            if sum(Gra_rrm==y)>0
                
                % select the output that is the closest to the accept pointer
                
                % temp_y returns all the outputs that send grants to input y
                temp_y = (Gra_rrm==y);
                
				acc_pointer(y)=0;
				
                % randomly select an output
                rnd_y = rand(1,N).*temp_y;   % generate random numbers      
                temp_o = find(rnd_y==max(rnd_y)); % pick the one with the largest random number
                
               
                % shift the queue forward by one packet
                % delete the transmitted packet
                temp_q = input_voq(y,temp_o,[2:end]);
                input_voq(y,temp_o,[1:end-1]) = temp_q;
                
                % decrease the queue length by 1
                voq_l(y,temp_o) = voq_l(y,temp_o) - 1;
                
                % calculate delay
                total_delay = total_delay + (t-voq_arr(y,temp_o,1));
                
                
                % delete the arrival time of the transmitted packet
                temp_p = voq_arr(y,temp_o,[2:end]);
                voq_arr(y,temp_o,[1:end-1]) = temp_p;
                
                % calculate total packets transmitted
                total_pkt = total_pkt + 1;
                
            end
        end
        
        t = t+1;        
    end
    
    % average cell latency
    Avg_delay_rrm(r) = total_delay/total_pkt;
    
end

out=[rho;Avg_delay_rrm];
% save rrm.mat out;
% 
% figure(1)
% semilogy(rho,Avg_delay_rrm);
% axis([0.2 1 1e-1 1e+3]);

