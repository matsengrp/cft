# Acquiring data from Immport

Go to [Immport's website](https://immport.niaid.nih.gov/immportWeb/home/home.do?loginType=full).
At the top, hover over "Access Data" and click "Download Data" from the drop-down menu.
To download data, you will need both an account and Aspera Connect installed.

For your account, your password cannot have any special characters for some reason.

For Aspera Connect, download a `.tar.gz` file---this contains a single `.sh` file that is incidentally the same size as the `.tar.gz` file.
Attempting to run as root will yield:

```
Installing Aspera Connect

This script cannot be run as root, Aspera Connect must be installed per user.
```

Fair enough.
Running without root gives:

```
Installing Aspera Connect

Deploying Aspera Connect (/home/dshaw/.aspera/connect) for the current user only.
Restart firefox manually to load the Aspera Connect plug-in

Install complete.
```

So in addition to a poor-entropy password, you also need Firefox, and you need to restart it.
I was not running Firefox at the time, so I opened it.
Aspera Connect was not working; it was only upon restarting firefox after opening it for the first time that Aspera Connect worked.

There is no way to point Aspera Connect to a specific directory path to place downloaded items, so your local machine must have at least the amount of free space as you are attempting to download.
