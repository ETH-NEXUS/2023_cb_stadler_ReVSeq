<script setup lang="ts"></script>

<template>
  <section class="tw-pt-12 tw-pb-14">
    <div class="tw-container tw-px-4 tw-mx-auto">
      <div class="tw-max-w-3xl tw-mx-auto">
        <p class="tw-text-lg tw-leading-8 tw-mb-6">
          This is the interface of the ReVSeq project that aims at providing an overview on the genetic
          diversity of respiratory viruses in Switzerland as well as fostering an understanding of their
          evolution. Samples from patients with general respiratory virus symptoms are collected at Viollier
          AG (as part of the routine diagnostics).
        </p>
        <img class="tw-block tw-w-full tw-rounded-md mb-2" src="../assets/research.jpg" alt="" />

        <p class="tw-text-lg tw-leading-8 tw-my-6">
          Viral RNA is extracted at
          <a class="tw-text-blue-600 hover:tw-text-purple-700" href="https://www.viollier.ch/de">
            Viollier AG
          </a>
          from the samples and undergoes high-throughput sequencing. The sequencing data is then analyzed at
          ETH Zürich by NEXUS
          <a class="tw-text-blue-600 hover:tw-text-purple-700" href="https://www.nexus.ethz.ch/">
            Personalized Health Technologies
          </a>
          and the
          <a class="tw-text-blue-600 hover:tw-text-purple-700" href="https://bsse.ethz.ch/cevo">
            Computational Evolution group
          </a>
          from Tanja Stadler.
        </p>
        <p class="tw-text-lg tw-leading-8 tw-my-6">
          First, all traces of human DNA sequence are removed from the samples. Second, the obtained sequence
          data is matched to reference sequences of known respiratory viruses, such as influenza or RSV.
          Third, virus consensus sequences are generated, and finally the type of virus causing the infection
          in each sample is identified. Virus sequence data and consensus sequences are shared with the
          public, by uploading data to
          <a class="tw-text-blue-600 hover:tw-text-purple-700" href="https://www.ebi.ac.uk/ena/browser/home">
            ENA
          </a>
        </p>
      </div>
    </div>
  </section>
</template>

<style scoped></style>
