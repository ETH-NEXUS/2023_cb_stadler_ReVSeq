<script setup lang="ts">
import {columnsSampleCount, columnsMetadata} from 'components/data'
import {createSampleCountsData, createMetadataRows} from 'components/helpers'
import {useCoreStore} from 'stores/core'

const coreStore = useCoreStore()
</script>

<template>
  <div class="q-pa-md" v-if="coreStore.sampleCounts.length > 0">
    <h3 class="text-h4 text-center q-pb-md">Counts</h3>
    <q-table
      :rows-per-page-options="[10, 20, 50, 100]"
      rows-per-page="20"
      flat
      bordered
      :rows="createSampleCountsData(coreStore.sampleCounts)"
      :columns="columnsSampleCount"
      row-key="name"></q-table>
  </div>
  <div class="q-pa-md" v-if="coreStore.metadata.length > 0">
    <h3 class="text-h4 text-center q-pb-md">Metadata</h3>
    <q-table
      :rows-per-page-options="[50, 100]"
      rows-per-page="50"
      flat
      bordered
      :rows="createMetadataRows(coreStore.metadata)"
      :columns="columnsMetadata"
      row-key="name"></q-table>
  </div>
</template>

<style scoped></style>
