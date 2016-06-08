import kddg.api.settings as settings
sys_settings = settings.load()

print(sys_settings.database)
print(sys_settings.database.username)
