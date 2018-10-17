from datetime import datetime
import os

class log_command(object):
    def __init__(self, command, from_function, th, function_class):
        self.command = command
        self.from_function = from_function
        self.th = th
        self.logs = {'function': "", 'command': "", 'start_time': "", 'end_time': "", 'threads': "", 'success': 0}
        self.f_class = function_class
        self.log = ""
        self.system_command_send()

    def system_command_send(self):
        self.log = self.f_class + ","
        self.log += self.from_function + ","
        self.log += self.th + ","
        self.log += str(datetime.now()) + ","

        try:
            os.system(self.command)
            self.log += str(datetime.now()) + ","
            self.log += "success,"
            self.log += self.command + "\n"
            self.write_logs(self.log)

        except:
            self.log += str(datetime.now()) + ","
            self.log += "failed,"
            self.log += self.command + "\n"
            self.write_logs(self.log)
            return self.from_function + " give error with this command -> " + self.command

    def write_logs(self, log):
        with open('log_file.txt', 'a') as file:
            file.write(log)
