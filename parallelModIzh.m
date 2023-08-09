configCluster("cm2_tiny", "cm2_tiny");
c = parcluster;
c.AdditionalProperties.PartitionName = 'serial_std';
c.AdditionalProperties.Cluster = 'serial_std';
c.AdditionalProperties.WallTime = '01:00:00';
c.AdditionalProperties.MemUsage = '20G';
c.AdditionalProperties.AdditionalSubmitArgs = '--cpus-per-task=4';

jobdir = fullfile(getenv('SCRATCH'), 'MdcsDataLocation/coolmuc', version('-release'));
if ~exist(jobdir), mkdir(jobdir); end
c.JobStorageLocation = jobdir;
c.saveProfile;

job = c.batch(@modIzhikevich, 0, {}, 'Pool', 2);

