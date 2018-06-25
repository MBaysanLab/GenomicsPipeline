from datetime import datetime
import os

class log_command(object):
    def __init__(self, command, from_function, th):
        self.command = command
        self.from_function = from_function
        self.th = th
        self.logs = {'function': "", 'command': "", 'start_time': "", 'end_time': "", 'threads': "", 'success': 0}
        self.system_command_send()

    def system_command_send(self):
        self.logs["function"] = self.from_function
        self.logs["command"] = self.command
        self.logs["start_time"] = str(datetime.now())
        self.logs["threads"] = self.th

        try:
            os.system(self.command)
            self.logs["end_time"] = str(datetime.now())
            self.logs["success"] = 1
            self.write_logs(self.logs)

        except:
            self.logs["end_time"] = str(datetime.now())
            self.logs["success"] = 0
            self.write_logs(self.logs)
            return self.from_function + " give error with this command -> " + self.command

    def write_logs(self, log):
        import json
        with open('log_file.txt', 'a') as file:
            file.write(json.dumps(log))
            file.write(",")