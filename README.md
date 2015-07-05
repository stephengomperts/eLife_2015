# eLife_2015

This repository contains matlab code accompanying the paper "VTA neurons coordinate with the hippocampal reactivation of spatial experience" by Stephen Gomperts, Fabian Kloosterman and Matthew Wilson.

Several utility functions that are used throughout the code can be found here: https://github.com/fkloosterman/matlab
Methods for neural decoding and replay trajectory fitting are described elsewhere (Davidson et al., 2009; Kloosterman, 2012; Kloosterman et al., 2014) and are not included here.

---

decoding.m

This file contains the skeleton code for showing how position dat is processed and how the decoders are set up for decoding.
For implementation of the decoding algorithm we refer to Kloosterman et al., 2014.

---

replay_detection_and_shuffles.m

This code evaluates SPW-R events to find significant replay events on the spatial working memory task. Specifically, 
it computes the probability that linear spatial trajectories of replayed position during SPW-R events are significantly different from chance. 
The linear trajectories of four trajectories are evaluated, given the four trajectories possible on the SWM task.
It takes as input an array of spatially decoded, concatenated candidate replay events, comprising a spatial pdf (10cm bins) x time (25msec bins), 
as well as the spatial bins of each possible trajectory on the SWM.
It uses a line finding algorithm based on Davidson et al. 2009 and further described in Kloosterman, 2012, that uses the radon
transform to find the best line in the spatial pdf of each candidate replay event.

---

rewardsitebias.m

This code is used to compute the reward site bias of replay events and of VTA units coordinating with replay events. 
It takes as input an array of spatially decoded, concatenated candidate replay events, comprising a spatial pdf (10cm bins) x time (25msec bins). 
It outputs the reward site bias for replay and for VTA spike-associated replay.

---

reward_site_bias_fwd_rvs.m

This code is used to compute the reward site bias of forward and reverse replayed spatial content, of centrifugal and centripetal replayed spatial content, and of VTA units coordinating with them.
It takes as input an array of spatially decoded, concatenated candidate replay events, comprising a spatial pdf (10cm bins) x time (25msec bins); 
a vector of replay slopes (m/s); an array of direction decoded, concatenated candidate replay events (outbound and inbound probability) x time (25ms bins), 
a vector of the rat's position at each replay event, and user-specified strings to assess (1) outbound/inbound decoded direction relative to track coordinates and (2) outgoing/incoming replay events relative to the animal's position. 
It outputs the reward site bias for replay and for VTA spike-associated replay, for each specific replayed spatial content 
(outbound & outgoing, outbound & incoming, inbound & outgoing, inbound & incoming), to derive the VTA unit associated reward site bias of forward, reverse, centrifugal and centripetal replayed spatial content.

---

theta_analysis.m

This code computes the hippocampal theta phase associated with each VTA spike
and the mean preferred phase and circular concentration of each VTA unit at hippocampal theta, 
used in Figure 6.
It uses the filter coefficients generated in filter_theta.m using the Matlab filter design tool.

---

sleep_scoring.m

This code identifies sleep segments and sleep-associated SPW-R events, used in Figure 7. 
It uses the filter coefficients generated in filter_delta.m and filter_theta.m using the Matlab filter design tool.
