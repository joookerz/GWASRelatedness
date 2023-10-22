print("This file was created within RStudio")

print("And now it lives on GitHub")

install.packages("usethis")
library(usethis)
git_sitrep()

usethis::create_github_token()
gitcreds::gitcreds_set()

use_git_config(user.name = "joookerz", user.email = "1435138480@qq.com")
#test3