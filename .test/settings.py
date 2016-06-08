import sys
sys.path.insert(0, '../..')

import ddg.ddglib.settings as settings
sys_settings = settings.load()

print(sys_settings.database)
print(sys_settings.database.username)
