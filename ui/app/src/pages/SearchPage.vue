<script setup lang="ts">
import SelectCard from 'components/SelectCard.vue'
import DataTables from 'components/DataTables.vue'
</script>

<template>
  <section>
    <div class="row justify-center tw-mb-16 flex-col">
      <div class="q-pa-md col-lg-4 col-md-6 col-sm-12 col-xs-12">
        <SelectCard />
      </div>

      <DataTables />
    </div>
  </section>
</template>
