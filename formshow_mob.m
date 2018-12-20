function html=formshow_mob(header,config)

data=header.Content;
data.seq=upper(data.seq(data.seq=='A'|data.seq=='T'|data.seq=='G'|data.seq=='C'|data.seq=='a'|data.seq=='t'|data.seq=='g'|data.seq=='c')); %data.seq=upper(regexprep(data.seq,'\r\n|\n|\r| ','')); %remove whitespace and end of line ncbi problem
data.W=str2double(header.Content.W);
data.T=str2double(header.Content.T);

%check if input data is ok and if not make warning html
if isempty(data.email) ||  ~(length(data.seq)==230) || isempty(data.origin) || isempty(data.W) || isempty(data.T)
    
    %create html
    html='<html><head></head><body>';
    html=[html 'Please provide all data correctly <br><br>'];
    html=[html '<br><br> <a href="form.html">Back to previous page</a><br>'];
    html=[html '<br><a href="../plasmid.html">Back to main page</a><br>'];
    html=[html '</body></html>'];

else
dir1 = pwd;
cd([dir1,'/www/MOB/formtest/'])
    
%save job number
load job.mat job
job=job+1;
save job.mat job

%save all data according to job number
name=sprintf('job%d_data.mat',job);
data.job=job;
save(name,'data');

%run computation
out=predictHostRange(job,data.seq,data.W,data.T);
data.out=out;
save(name,'data');
% name=sprintf('job%d_out.csv',job);
% csvwrite(name,out)

%publish or email results
% subject=sprintf('TIDD prediction server job %d results',job);
% message='Hello, your job has completed and the results are attached. Thank you for using our server, yours sincerely, Lapanjes Lab team';
% attachments{1}=name;
% attachments{2}='plot.jpg';
% 
% user_name = 'aa';                          % My Username
% password = 'aa';      % My E-Mail Password
% smtp_server = 'smtp.gmail.com';               % My SMTP Server
% setpref('Internet', 'SMTP_Username', user_name);
% setpref('Internet', 'SMTP_Password', password);
% setpref('Internet', 'SMTP_Server', smtp_server);
% props = java.lang.System.getProperties;
% props.setProperty('mail.smtp.auth', 'true');
% props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
% props.setProperty('mail.smtp.socketFactory.port', '465');
% 
% sendmail(email,subject,message,attachments)

cd(dir1)

%create html
str1 = ['<br><br>MOB group: ',out.mob,'<br>'];
str2 = ['<br><br>Repertoire of potential hosts: ',out.host,'<br>']; % display results

html='<html><head></head><body>';
html=[html 'The results for your query are: <br><br>'];
html=[html str1]; % display results
html=[html str2]; % display results
html=[html '<br><br>Legend for results of host range*: <br>'];
html=[html '<img src="Host_range.png" width="28%"/>'];
html=[html '<br><br>* Taxonomic families of Gammaproteobacteria where either the reference plasmid or other pasmids of the same Inc/REP group have been found.<br>'];
html=[html '<br><br> <a href="form.html">Back to previous page</a><br>'];
html=[html '<br><a href="../plasmid.html">Back to main page</a><br>'];
html=[html '</body></html>'];

% html='<html><head></head><body>';
% html=[html 'This html is produced by formshow.m <br><br>'];
% html=[html 'Your user Credentials<br>'];
% ContentForm=headers.Content;
% fn=fieldnames(ContentForm);
% for i=1:length(fn)
% 	html=[html 'Your ' fn{i} ' is ' ContentForm.(fn{i}) '<br>'];
% end
% html=[html '<br><br> <a href="form.html">Back to Previous Page</a><br>'];
% html=[html '</body></html>'];

end 
