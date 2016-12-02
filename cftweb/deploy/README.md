# Overview

The CFT dynamic web service needs to be hosted on some machine within
the Hutch.  There is a open-source platform called `proxmox` that is
used within the Hutch to create virtual machines on-demand, without
having to go through a help ticket.

You can read more about how to create a VM inside the Hutch on the SciComp Wiki.
Here is one link that I found particularly helpful:

https://teams.fhcrc.org/sites/citwiki/SciComp/Pages/How%20to%20run%20your%20own%20Docker%20host.aspx

## Creating the cftweb server

The prox commands to create the cftweb server are in `delpoy.sh`.  The
script will create a server called `cftweb4.fhcrc.org`.  The server
will be accessible only inside the Hutch, not from the internet.

The cftweb server is run as a docker container on the cftweb VM.

## `prox` commands that are useful to know

A proxmox dashboard sowing the status of all VMs is available at
https://proxa1.fhcrc.org:8006/ Log in with your hutch net id (select
Active Directory authentication)

In addition you can interact with the proxmox server through
commandline tools that are available on all the gizmo/scicomp machines.

   - prox list - see a list of your VMs
   - prox ssh <machine> - ssh in to the VM named <machine>
   - prox stop <machine> - shut down the VM named <machine>
   - prox start <machine> - boot the VM named <machine>
   - prox destroy <machine> - *READ CAVEAT BELOW*  destroy the VM named <machine>


## IMPORTANT - shutting down virtual machines.

When creating a new VM with `prox`, a DNS entry will be
autmatically created.  The DNS name is bound to a DHCP lease which 3
days right now.  If you destroy a VM, the DNS entry still hang around
pointing to the previously assigned ip address.  The DNS entry will
expire after three days, but in the meantime you cannot reuse that
name for another VM.

There are 3 ways for you to reuse the name of a VM:

1. Use shutdown –h inside the VM and then destroy it with `prox destroy`
2. Send a ticket to helpdesk and ask to delete the DNS name
3. You can pick a new hostname temporarily, then wait 3 days, rename
   it to the old name in the hosts file and reboot. It should pick up
   the old name.

